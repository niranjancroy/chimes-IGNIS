# run as >>> execfile('bhfeedback_newsnaplist.py')

import numpy as np
import gadget as g
import h5py as h5py
from glob import glob
import os as os
from daa_constants import *
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM, z_at_value
np.set_printoptions(linewidth=180)



# --- original snapshot
sim_tag = 'h113_HR'
#snapnum = 153
snapnum = 214
#snapnum = 248


sim_name = sim_tag + '_sn' + g.snap_ext(snapnum)
inidir = '/simons/scratch/dangles/FIRE/bhfeedback/' + sim_name + '/ICs'
outdir = inidir


# --- new ICs
#dt_Myr = np.array( [ 0.01, 0.1, 1. ] )
#t_int = np.array( [ 0., 1., 20., 121. ] )

dt_Myr = np.array( [ 0.01 ] )
t_int = np.array( [ 0., 10. ] )


# --- create output times file:
P = g.readsnap(inidir,0,5,cosmological=1)
z_begin = P['redshift']
cosmo = FlatLambdaCDM(H0=P['hubble']*100, Om0=P['omega_matter'])
t_begin = cosmo.age(z_begin).to(u.Myr)


t_all = np.array([])
for i in range(dt_Myr.size):    
   t_all = np.append( t_all, np.arange( t_int[i], t_int[i+1], dt_Myr[i] ) )

dt_all = t_all[1:] - t_all[:-1]

np.set_printoptions(precision=2, suppress=True)
print '\nt_all =', t_all
print '\ndt_all =', dt_all
print '\nnsnap =', t_all.size

times_all = t_all * u.Myr + t_begin

z_all = np.array([])
for i in range(times_all.size):
   z_all = np.append( z_all, z_at_value(cosmo.age,times_all[i]) )
a_all = 1./(1.+z_all)

if dt_Myr.size > 1:
  filename = outdir+'/snaplist_sn%s_z%.2f_%.2f-%.0fMyr.txt' % ( snapnum, z_begin, dt_Myr[0],dt_Myr[-1] )
else:
  filename = outdir+'/snaplist_sn%s_z%.2f_%.3fMyr.txt' % ( snapnum, z_begin, dt_Myr[0] )

#np.savetxt( filename, a_all, fmt='%10.8f' )
print '\n',filename


print '\nDone!!'


