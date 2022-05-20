from __future__ import division,print_function

import numpy as np

def nfw(log_r, log_norm, log_rscale):
	xc = log_r - log_rscale;
	return log_norm - xc - 2.*np.log10(1.+10.**xc);
	
def double_pwrlaw(log_x, log_lstar, log_phistar, faint_slope, bright_slope):
    x = 10.**(log_x-log_lstar);
    return log_phistar - np.log10(x**faint_slope + x**bright_slope);

def beta_of_n_sersic(n):
    # beta is beta**1/n, which is much better behaved: get excellent 
    #   fit to the full numerical solution (to <1% for 0.1<n_s<20) with:
    d0=0.163382; b0=1.9994830; c0=-0.32670002; n0=4.;
    b=3.*n*n*np.exp(-(n/(1.5*d0))**n0) + (c0+b0*n)/(1.+10.*np.exp(-(n/d0)**n0)/((n/0.1)**n0))
    return b;

def schechter(x, alpha, L_star, phi_star):
    # returns the linear value of Schechter function in luminosity-space
    # (but note, STILL returning dphi/dlogL, *NOT* dphi/dL)
    return np.log(10.)*phi_star * ((x/L_star)**(1.-alpha)) * EXP(-1.0*x/L_star);
    
def sersic(log_r, log_I_e, log_r_e, n):
    r = 10.**log_r
    r_e = 10.**log_r_e;
    b_n = beta_of_n_sersic(n);
    return log_I_e - (b_n/np.log(10.))*( (r/r_e)**(1./n) - 1.);
    
def sersic_disk(log_r, log_I_e, log_r_e, n, log_I_e_2, log_r_e_2, n_2):
    r = 10.**log_r

    r_e = 10.**log_r_e;
    b_n = beta_of_n_sersic(n);
    profile_1 = log_I_e - (b_n/np.log(10.))*( (r/r_e)**(1./n) - 1.);

    r_e = 10.**log_r_e_2;
    b_n = beta_of_n_sersic(n_2);
    profile_2 = log_I_e_2 - (b_n/np.log(10.))*( (r/r_e)**(1./n) - 1.);

    return np.log10( 10.**profile_1 + 10.**profile_2 );