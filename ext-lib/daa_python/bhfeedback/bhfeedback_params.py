# >>> execfile('bhfeedback_params.py')

import numpy as np
import gadget as g
import daa_lib as daa
from daa_constants import *


c = SPEED_OF_LIGHT / CM_PER_KM   # speed of light in km/s   
rad_eff = 0.1


# --- BH parameters
print '\n---- BH parameters:'
Mbh = 1e8                             # Msun
edd_ratio = 1.
Mdot = daa.MdotEddington( Mbh )       # Msun/year
Lbol = rad_eff * (Mdot*MSUN/SEC_PER_YEAR) * SPEED_OF_LIGHT**2          # erg/s

print 'Mbh = {:.2g} Msun,  edd_ratio = {:.2f},  Mdot = {:.3g} Msun/yr,  Lbol = {:.3g} erg/s'.format( Mbh, edd_ratio, Mdot, Lbol )

# --- BAL parameters (inner accretion disk winds)
print '\n---- BAL wind parameters:'
#for BAL_v, BAL_facc in zip([3e4,3e4,3e4],[1./11.,1./3.,0.5]):
for BAL_v, BAL_facc in zip([3e4,3e4,3e4,3e4,3e4],[0.50000000,0.33333333,0.20000000,0.11111111,0.09090909]):
   #BAL_v = 3e4     # km/s
   #BAL_facc = 0.5 
   BAL_MdotOut_over_Mdot = (1.-BAL_facc)/BAL_facc
   BAL_MdotOut = BAL_MdotOut_over_Mdot * Mdot
   BAL_Edot = 0.5 * (BAL_MdotOut*MSUN/SEC_PER_YEAR) * (BAL_v*CM_PER_KM)**2
   BAL_Edot_over_Lbol = BAL_Edot / Lbol
   BAL_MomFlux_over_Lbolflux = (BAL_MdotOut*MSUN/SEC_PER_YEAR) * (BAL_v*CM_PER_KM) / (Lbol/SPEED_OF_LIGHT)

   print 'BAL_v = {:.2g} km/s,  BAL_facc = {:.3f},  BAL_MdotOut_over_Mdot = {:.3f},  BAL_MdotOut = {:.3g} Msun/yr,  BAL_Edot = {:.3g} erg/s,  BAL_Edot_over_Lbol = {:.3f},  BAL_MomFlux_over_Lbolflux = {:.2f}'.format( BAL_v, BAL_facc, BAL_MdotOut_over_Mdot, BAL_MdotOut, BAL_Edot, BAL_Edot_over_Lbol, BAL_MomFlux_over_Lbolflux )

#print '\n'
#for BAL_MdotOut_over_Mdot in [1., 2., 4., 6., 8., 10.]:
#   BAL_facc = 1. / (1. + BAL_MdotOut_over_Mdot)
#   print 'BAL_MdotOut_over_Mdot = {:.3f}, BAL_facc = {:.8f}'.format(BAL_MdotOut_over_Mdot, BAL_facc)

# --- particle SPAWNING parameters
print '\n---- particle SPAWNING parameters:'
SPAWN_v_shock = BAL_v / 4.
SPAWN_T_shock = 1.2e10 * ( BAL_v/(0.1*c) )**2

print 'SPAWN_v_shock = {:.2g} km/s,  SPAWN_T_shock = {:.3g} K'.format( SPAWN_v_shock, SPAWN_T_shock, )


# --- KICK wind parameters
print '\n---- KICK wind parameters:'
KICK_Edot = BAL_Edot   # same kinetic energy as BAL winds
KICK_Edot_over_Lbol = KICK_Edot / Lbol
for KICK_MomFlux_over_Lbolflux in [ 1., 2., 5., 10., 20. ]:     # fix the effective momentum flux
   KICK_MomFlux = KICK_MomFlux_over_Lbolflux * (Lbol/SPEED_OF_LIGHT)
   KICK_v = 2. * KICK_Edot / KICK_MomFlux / CM_PER_KM                          # km/s
   KICK_MdotOut = KICK_MomFlux**2. / (2.*KICK_Edot) * SEC_PER_YEAR / MSUN      # Msun/year
   KICK_MdotOut_over_Mdot = KICK_MdotOut / Mdot
   print 'KICK_Edot_over_Lbol = {:.2f},  KICK_MomFlux_over_Lbolflux = {:.2f},  KICK_v = {:.4g} km/s,  KICK_MdotOut_over_Mdot = {:.2f},  KICK_MdotOut = {:.3g} Msun/yr'.format(
          KICK_Edot_over_Lbol, KICK_MomFlux_over_Lbolflux, KICK_v, KICK_MdotOut_over_Mdot, KICK_MdotOut )

# --- SIMBA wind/jet parameters
print '\n---- SIMBA wind/jet parameters:'
for SIMBA_v in [500., 1000., 7000., 10000.]:     # fix the wind velocity in km/s
   SIMBA_MomFlux_over_Lbolflux = 20.
   SIMBA_MomFlux = SIMBA_MomFlux_over_Lbolflux * (Lbol/SPEED_OF_LIGHT)              # in cgs
   SIMBA_MdotOut = SIMBA_MomFlux / (SIMBA_v * CM_PER_KM)                            # in cgs
   SIMBA_MdotOut_over_Mdot = SIMBA_MdotOut * (SEC_PER_YEAR/MSUN) / Mdot
   SIMBA_Edot = 0.5 * SIMBA_MdotOut * (SIMBA_v*CM_PER_KM)**2                        # in cgs
   SIMBA_Edot_over_Lbol = SIMBA_Edot / Lbol
   print 'SIMBA_Edot_over_Lbol = {:.2f},  SIMBA_MomFlux_over_Lbolflux = {:.2f},  SIMBA_v = {:.4g} km/s,  SIMBA_MdotOut_over_Mdot = {:.2f},  SIMBA_MdotOut = {:.3g} Msun/yr'.format(
          SIMBA_Edot_over_Lbol, SIMBA_MomFlux_over_Lbolflux, SIMBA_v, SIMBA_MdotOut_over_Mdot, SIMBA_MdotOut )


# --- time/length scales
v_test = np.array([1e2,1e3,1e4])
print '\n---- time/length scales:'
print 'vel = ', v_test, ' km/s   x  1 Myr  = ', v_test * 1e6 * SEC_PER_YEAR * CM_PER_KM / CM_PER_KPC, ' kpc'

# --- Schwarzschild radius:

Rs = 2. * G_UNIV * Mbh * MSUN / SPEED_OF_LIGHT**2 / (CM_PER_KPC/1e3)
print '\n---- Schwarzschild radius::'
print 'Rs =', Rs, ' pc'




"""
mom_flux_ALL = [  10, 10, 10, 10 ]
v_kpc_ALL    = [ 5e2, 1e3, 5e3, 1e4 ]
e_m_ALL      = [ 0.1, 0.1, 0.1, 0.1 ]


for i in range(len(mom_flux_ALL)):

   mom_flux = mom_flux_ALL[i] 
   v_kpc = v_kpc_ALL[i]
   e_m = e_m_ALL[i]

   # --- derived quentities ---

   M_kpc_over_Mdot = mom_flux * rad_eff * ( c / v_kpc )

   f_acc = 1. / ( 1. + M_kpc_over_Mdot )

   e_kin = (1./(2.*rad_eff)) * ( v_kpc/c )**2 * ( 1. - f_acc ) / f_acc
   e_kin_tmp = 0.5 * mom_flux * v_kpc / c

#   M_0_over_Mdot = ( 1. - e_m ) / e_m
#   v_0 = np.sqrt( M_kpc_over_Mdot / M_0_over_Mdot ) * v_kpc               # energy conservation

#   print '\n-------------------------------------\ninput parameters:'
#   print 'v_kpc =', v_kpc, '   e_m = ', e_m, '   f_acc =', f_acc
#   print '\ncorrespond to...'
#   print 'mom_flux =', mom_flux, '   v_0 =', v_0, '   M_kpc/Mdot = ', M_kpc_over_Mdot, '    e_kin = ', e_kin, '\n-------------------------------------\n'


   print ' \nINPUT:  mom_flux = Pdot/(Lbol/c) =', mom_flux, '  v_kpc =', v_kpc, 'km/s'
   print 'OUTPUT:  M_out/Mdot =', M_kpc_over_Mdot, '  e_kin = Edot/Lbol = ', e_kin, '  f_acc =', f_acc

"""

print "\n\nDone!\n"



