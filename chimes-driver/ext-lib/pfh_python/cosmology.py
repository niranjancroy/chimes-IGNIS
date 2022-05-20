from __future__ import division,print_function

import numpy as np
import scipy.misc
from astropy.cosmology import WMAP9 as cosmo
""" 
astropy.cosmology sets dictionary with :
    ==========  =====================================
    Oc0         Omega cold dark matter at z=0
    Ob0         Omega baryon at z=0
    Om0         Omega matter at z=0
    flat        Is this assumed flat?  If not, Ode0 must be specifiec
    Ode0        Omega dark energy at z=0 if flat is False
    H0          Hubble parameter at z=0 in km/s/Mpc
    n           Density perturbation spectral index
    Tcmb0       Current temperature of the CMB
    Neff        Effective number of neutrino species
    sigma8      Density perturbation amplitude
    tau         Ionisation optical depth
    z_reion     Redshift of hydrogen reionisation
    t0          Age of the universe in Gyr
    reference   Reference for the parameters
    ==========  =====================================
"""


## returns universal baryon fraction
def baryon_fraction():
    return cosmo.Ob0/cosmo.Om0


## compute the critical overdensity for collapse vs. redshift z
def collapse_threshold( z ):
    delta_c_0 = 0.15*((12.*np.pi)**(2./3.))*((cosmo.Om(z))**(0.0055))     
    return delta_c_0/growth_factor(z)


## return the dimensionless time element per redshift element
def dt_dz( z ):
    return 1./((1.+z) * cosmo.efunc(z))


## returns the mass density (not critical density) of the universe
##   if 'cgs' is set, units are cgs, otherwise M_sun/Mpc^3
##   (physical, no h factors here)
def mass_density( z , cgs=0):
    rho_c = cosmo.critical_density(z).value
    rho_m = cosmo.Om(z) * rho_c
    if(cgs==1):
        return rho_m
    else:
        return rho_m / (1.989e33/((3.086e24)**3.))


## the linear growth factor D(z), normalized to D(z=0)=1, evaluated at a redshift z
def growth_factor( z , exact=False ):
    if (exact==True):
        return -1; ## not imported yet, from pfh idl library
    ## use Carroll et al. 1992 approximation (good to ~1%):
    omz = cosmo.Om(z)
    olz = cosmo.Ode(z)
    gz  = (5./2.)*omz/(omz**(4./7.) - olz + (1.+omz/2.)*(1.+olz/70.))
    omz = cosmo.Om0
    olz = cosmo.Ode0
    g0  = (5./2.)*omz/(omz**(4./7.) - olz + (1.+omz/2.)*(1.+olz/70.))
    Dz  = (gz/g0)/(1.+z)
    return Dz


## returns halo r_vir (*PHYSICAL*) as a function of mass & redshift
##   (r_vir in h^-1 kpc for m_halo in h^-1 M_sun)
##   (bryan & norman 1998 by way of barkana & loeb 2001)
def rvir_of_mhalo( log10_M_halo, z ):
    omz = cosmo.Om(z)
    d = omz-1.
    dc = 18.*np.pi*np.pi+82.*d-39.*d*d
    f = (cosmo.Om0 / omz) * (dc/(18.*np.pi*np.pi))
    rv = 0.784*10**((log10_M_halo-8.)/3.) * (f**(-1./3.)) * (((1.+z)/10.)**(-1.))
    return rv


## returns halo v_c as a function of mass & redshift
##   (vc in km/s for m_halo in h^-1 M_sun)
##   (bryan & norman 1998 by way of barkana & loeb 2001)
def vc_of_mhalo( log10_M_halo, z ):
    omz = cosmo.Om(z)
    d = omz-1.
    dc = 18.*np.pi*np.pi+82.*d-39.*d*d
    f = (cosmo.Om0 / omz) * (dc/(18.*np.pi*np.pi))
    vc = 23.4*10**((log10_M_halo-8.)/3.) * (f**(1./6.)) * (((1.+z)/10.)**(1./2.))
    return vc


## returns halo virial temperature as a function of mass & redshift
##   (t_vir in K for m_halo in h^-1 M_sun)
##   (bryan & norman 1998 by way of barkana & loeb 2001)
def tvir_of_mhalo( log10_M_halo, z ):
    omz = cosmo.Om(z)
    d = omz-1.
    dc = 18.*np.pi*np.pi+82.*d-39.*d*d
    f = (cosmo.Om0 / omz) * (dc/(18.*np.pi*np.pi))
    mu = 0.6 ## mean molecular weight -- canonical number for pristine gas
    return 1.98e4 * (mu/0.6) * 10.**((log10_M_halo-8.)*2./3.) * (f**(1./3.)) * ((1.+z)/10.)


## returns halo concentration parameter as a function of mass & redshift
##    (halo mass in units of h^-1 M_sun)
##  - could replace with more detailed (bullock et al. 2001) model, but fitting fn 
##      below from prada et al. (and others) works just as well for our purposes and 
##      extrapolates more stably to very high masses (less so very low masses?)
def concentration_of_mhalo( log10_M_halo, z ):
    a = 1./(1.+z)
    x = (cosmo.Ode0/cosmo.Om0)**(1./3.) * a
    B0=0.8767405803323409-0.08659849959728688*np.arctan(2.945952-6.948*x)
    B1=0.8329911404974490-0.11795353862546273*np.arctan(3.885036-7.386*x)
    sprime = B1 * growth_factor(z) * sigma_of_mhalo_z0(log10_M_halo)
    c_sprime = 2.881 * ((sprime/1.257)**(1.022) + 1.) * np.exp(0.060/(sprime*sprime))
    c = B0 * c_sprime
    return c    


## returns halo scale radius r_s (*PHYSICAL*) as a function of mass & redshift
##   (r_scale in h^-1 kpc for m_halo in h^-1 M_sun)
##    NFW-like profile assumed with r_scale = r_vir / concentration
def rscale_of_mhalo( log10_M_halo, z ):
    return rvir_of_mhalo(log10_M_halo,z)/concentration_of_mhalo(log10_M_halo,z)


## return the maximum circular velocity of a halo as a function of mass & redshift
##   (vc in km/s for m_halo in h^-1 M_sun)
def vmax_of_mhalo( log10_M_halo, z ):
    vv = vc_of_mhalo( log10_M_halo, z )
    cv = concentration_of_mhalo( log10_M_halo, z)
    Av = np.log(1.+cv) - cv/(1.+cv)
    vmax = vv * np.sqrt( 0.216 * cv / Av )
    return vmax


## function which returns mstar of the DM mass function 
##   (where sigma(M) = d_coll at the given redshift)
##   returns (log(Mstar)) in h^-1 M_sun
def mstar_dm( z ):
    d_c = collapse_threshold(z)
    m_grid = np.arange(2.,17.,0.05)
    s_grid = sigma_of_mhalo_z0(m_grid)
    m_star = np.interp(d_c,s_grid,m_grid)
    return m_star


## function to return the local sigma(M) as a function of Mhalo (in Msun/h)
##   !! designed to be done in h_inverse coordinates !! (but can switch to h's below)
def sigma_of_mhalo_z0( log10_M_halo , exact=False):
    if (exact==True):
        return -1; ## need to import still (in pfh idl libraries)
    ##
    ## here's the default -- the fitting formula from van den Bosch 2002 (mnras, 331, 98), 
    ##   modified with a tiny tweak to accommadate a non-unity spectral index, calibrated 
    ##   against both the WMAP1, WMAP3, and concordance cosmologies
    ##   -- gives s(M) to within a factor 0.005 error, for 10^6 < M < 10^16 
    ##      (overpredicts power at M > 10^16, underpredicts at < 10^6; although still 
    ##        good at ~5% level down to M ~ 10^2 )
    ##
    omega_matter = cosmo.Om0
    hubble = cosmo.H0.value/100.
    omega_baryon = cosmo.Ob0
    sigma_8 = 0.8
    spectral_index = 0.98
    gamma = omega_matter*hubble*np.exp(-omega_baryon * (1.+np.sqrt(2.*hubble)/omega_matter))
    c = 3.804e-4 
    x = (c*gamma*((10**(log10_M_halo/3.))))/((omega_matter)**(1./3.))
    g1 = 64.087 * ((1. + 1.074*(x**(0.3)) - 1.581*(x**(0.4)) + 0.954*(x**(0.5)) - 0.185*(x**(0.6)))**(-10.))
    x = (32.*gamma)
    g2 = 64.087 * ((1. + 1.074*(x**(0.3)) - 1.581*(x**(0.4)) + 0.954*(x**(0.5)) - 0.185*(x**(0.6)))**(-10.))
    f = (g1*g1)/(g2*g2)
    s = np.sqrt(f * sigma_8 * sigma_8)
    return s*(10**((log10_M_halo-14.09)*(1.00-spectral_index)/9.2))/(1.+(1.00-spectral_index)/9.2)


## function to return the bias of a halo of mass Mhalo (in Msun/h) at redshift z
def bias_of_mhalo( log10_M_halo, z ):
    """
    ## return analytical Sheth-Tormen bias for a given halo mass at some redshift
    a = 0.75
    p = 0.3
    delta_c = collapse_threshold(REDSHIFT)
    sigma_M = sigma_of_mhalo_z0(log_M_halo) * growth_factor(REDSHIFT)
    nu = delta_c/sigma_M
    b = 1. + (a*nu*nu - 1.)/delta_c + (2.*p)/(delta_c*(1. + (a*nu*nu)^p))
    """
    ## version from Croom et al. clustering paper -- Sheth, Mo, & Tormen (2001)
    a,b,c,sqrt_a = 0.707,0.5,0.6,0.840832920383
    delta_c = collapse_threshold(z)
    nu = delta_c / sigma_of_mhalo_z0(log10_M_halo)
    denom_1 = sqrt_a * delta_c*growth_factor(z)
    numer_1 = a*sqrt_a * nu*nu
    q = (a*nu*nu)**c
    numer_2 = numer_1 * b/q
    numer_3 = -q
    denom_3 = q + b*(1.-c)*(1.-c/2.)
    bias = 1. + (numer_1 + numer_2 + numer_3 / denom_3 ) / denom_1
    return bias 


## returns the average mass of the "main" progenitor halo 
##   to a z=0 halo of the given mass, at a given redshift
##  - calculated following the analytical solutions in 
##     Neistein et al. 2006
##
## (assumes all halo masses in h^-1 M_sun)
##
def descendant_mass( log10_m_halo, z, z_final=0. ):
    ## interpolate over a grid of possible masses
    mh_z0 = np.arange(3.,20.,0.25)
    progen_mass = progenitor_mass( mh_z0, z, z_init=0.0 )
    desc_mass_z0 = np.interp(log10_m_halo, progen_mass, mh_z0)
    desc_mass = desc_mass_z0
    if (z_final > 0.):
        if (z<=z_final):
            print('z must be > z_final!')
            return 0
        desc_mass = progenitor_mass( desc_mass_z0, z_final, z_init=0.0 )
    return desc_mass


def progenitor_mass( log10_m_halo, z , z_init=0.):
    log10_m_halo_z0 = log10_m_halo
    if (z_init > 0.):
        if (z<=z_init):
            print('z must be > z_init!')
            return 0
        log10_m_halo_z0 = descendant_mass( log10_m_halo, z_init )
    omega_matter = cosmo.Om0
    hubble = cosmo.H0.value/100.
    omega_baryon = cosmo.Ob0
    sigma_8 = 0.8
    gamma = omega_matter*hubble*np.exp(-omega_baryon * (1.+np.sqrt(2.*hubble)/omega_matter))
    w  = collapse_threshold(z) 
    w0 = collapse_threshold(0.)
    dw  = w - w0
    x   = 32. * gamma
    g   = 64.087 * ((1. + 1.074*(x**(0.3)) - 1.581*(x**(0.4)) + 0.954*(x**(0.5)) - 0.185*(x**(0.6)))**(-10.))
    f_1 = (sigma_8/g)
    x   = (gamma**3.) * (10.**(log10_m_halo_z0)) / omega_matter
    u   = x
    lnu = np.log(u)
    f_2 = (-6.92e-5)*((lnu)**(4.0)) + (5.0e-3)*((lnu)**(3.0)) + \
               (8.64e-2)*((lnu)**(2.0)) + (-12.66)*((lnu)) + 110.8
    f_3 = (1./f_1) * dw
    x   = f_3 + f_2 
    lnu = np.arange(40.,-10.,-0.1)*np.log(10.)
    f   = (-6.92e-5)*((lnu)**(4.0)) + (5.0e-3)*((lnu)**(3.0)) + \
               (8.64e-2)*((lnu)**(2.0)) + (-12.66)*((lnu)) + 110.8
    finv = np.exp(np.interp(x,f,lnu)) 
    ## remember, np.interp input 'x-axis' vector MUST be increasing!!
    mm  = omega_matter/(gamma**(3.)) * finv
    return np.log10(mm)


## returns the Sheth-Tormen
##   halo mass function ( dn/dlog(Mhalo) in Mpc^-3 log(Mhalo)^-1 -- no h's!) 
##   for a given mass (in h^-1 M_sun) and redshift
##   -- this is in comoving Mpc^-3, as standard for MFs
def mass_function( log10_M_halo, z, press_schechter=0 ):
    ## sheth-tormen MF parameters
    a,b,p = 0.707,0.3222,0.3
    rho_m = mass_density(0.) ## in comoving units, so total mass density is constant
    d_c = collapse_threshold(z)
    sig = sigma_of_mhalo_z0(log10_M_halo)
    hubble = cosmo.h
    dlns_dlogM = -scipy.misc.derivative(sigma_of_mhalo_z0,log10_M_halo)/sig
    x = np.sqrt(a) * d_c/sig
    f_x = b * np.sqrt(2./np.pi) * (1.+x**(-2.*p)) * x * np.exp(-x*x/2.)
    dn_dlogM = rho_m / ((10.**log10_M_halo)/hubble) * f_x * dlns_dlogM
    return dn_dlogM


## return age of universe (for a flat universe) to a given redshift
def age_of_universe(z,h=0.71,Omega_M=0.27):
    ## exact solution for a flat universe
    a=1./(1.+z); x=Omega_M/(1.-Omega_M) / (a*a*a);
    t=(2./(3.*np.sqrt(1.-Omega_M))) * np.log( np.sqrt(x) / (-1. + np.sqrt(1.+x)) );
    t *= 13.777 * (0.71/h); ## in Gyr
    return t;


## return lookback time (for a flat universe) to a given redshift
def lookback_time(z,h=0.71,Omega_M=0.27):
    return age_of_universe(0.,h=h,Omega_M=Omega_M)-age_of_universe(z,h=h,Omega_M=Omega_M);
