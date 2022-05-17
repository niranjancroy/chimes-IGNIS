# run as >>> execfile('bhfeedback_snapICs.py')

import numpy as np
import gadget as g
import h5py as h5py
from glob import glob
import os as os
from daa_constants import *
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM, z_at_value
np.set_printoptions(linewidth=150)


"""
fire_dir = '/projects/b1026/anglesd/FIRE/'
# --- original snapshot
inidir = fire_dir + 'h113_HR_sn1dy300ro100ss'
sim_tag = 'h113_HR'
snapnum = 152
# --- new ICs
logMbh = 8
sim_name = sim_tag + '_sn' + g.snap_ext(snapnum) 
outdir = fire_dir + 'bhfeedback/' + sim_name
"""


fire_dir = '/simons/scratch/dangles/SMAUG/h113_HR_sn152res251/ICs/'
# --- original snapshot
inidir = fire_dir + 's8'
snapnum = 0
# --- new ICs
#logMbh = 9
#sim_tag = 's9'
#logMbh = 7
#sim_tag = 's7'
logMbh = 6
sim_tag = 's6'


"""
fire_dir = '/simons/scratch/dangles/SMAUG/h113_HR_sn262res300/ICs/'
# --- original snapshot
inidir = fire_dir + 's8'
snapnum = 0
# --- new ICs
#logMbh = 9
#sim_tag = 's9'
#logMbh = 7
#sim_tag = 's7'
logMbh = 6
sim_tag = 's6'
"""


outdir = fire_dir + sim_tag + '/snapdir_000'

if not os.path.exists(outdir):
  os.makedirs(outdir)


P = g.readsnap(inidir,snapnum,5,cosmological=1)
ind_bh = np.argmax(P['mbh'])
bh_ID = P['id'][ind_bh]
z_begin = P['redshift']
a_begin = 1. / (1. + z_begin)

print '\nBH mass = ', P['mbh'][ind_bh] * UnitMass_in_Msun, '  redshift =', z_begin
print np.sort(P['mbh'])[::-1][:10]/P['mbh'][ind_bh], '\n'


snapfiles = glob( inidir + '/snapdir_' + g.snap_ext(snapnum) + '/snapshot*.hdf5')
for i, thisfile in enumerate(snapfiles):

   f = h5py.File( thisfile, 'r' )

   #fname = '/ics_' + sim_name + '_z%.7f_a%.7f' % ( z_begin, a_begin)
   #outfile = outdir + fname + '.' + thisfile.split('.')[-2] + '.hdf5'

   outfile = outdir + '/snapshot_000.' + thisfile.split('.')[-2] + '.hdf5'   

   print i, thisfile, outfile
   #continue


   ics = h5py.File(outfile, 'w')

   for key in f:
      f.copy(key, ics)
   for attr_key in f.attrs:
      ics.attrs[attr_key] = f.attrs[attr_key]

   # --- edit ics, we only keep the most massive BH
   del ics['Header'].attrs['NumPart_Total']
   del ics['Header'].attrs['NumPart_ThisFile']
   NumPart_Total = f['Header'].attrs['NumPart_Total'][:];  NumPart_Total[5] = 1
   NumPart_ThisFile = f['Header'].attrs['NumPart_ThisFile'][:]
   ics['Header'].attrs['NumPart_Total'] = NumPart_Total

   if 'PartType5' in f.keys():

     """
     ind_theBH = np.where( f['PartType5']['ParticleIDs'][:] == bh_ID )[0]
     if ind_theBH.size == 0:
       print '  no BHs'
       NumPart_ThisFile[5] = 0
       ics['Header'].attrs['NumPart_ThisFile'] = NumPart_ThisFile
       if 'PartType5' in ics.keys():
         del ics['PartType5']
     else:
       ind_theBH = ind_theBH[0]
       print '  BH_ID:  ', bh_ID, f['PartType5']['ParticleIDs'][ind_theBH], ics['PartType5']['ParticleIDs'][ind_theBH]
       NumPart_ThisFile[5] = 1
       ics['Header'].attrs['NumPart_ThisFile'] = NumPart_ThisFile
       for keyname in f['PartType5'].keys():           
          del ics['PartType5'][keyname]
          if len(f['PartType5'][keyname].shape) == 1:
            if keyname in ['Masses','BH_Mass']:
              ics['PartType5'][keyname] = 10.**logMbh / (UnitMass_in_Msun/ics['Header'].attrs['HubbleParam'])
            elif keyname in ['BH_Mass_AlphaDisk','BH_Mdot']:
              ics['PartType5'][keyname] = 0.
            else:
              ics['PartType5'][keyname] = f['PartType5'][keyname][ind_theBH]
          elif len(f['PartType5'][keyname].shape) == 2:
            ics['PartType5'][keyname] = f['PartType5'][keyname][ind_theBH,:]
          else:
            ics['PartType5'][keyname] = f['PartType5'][keyname]
     """
  
     print '  BH_ID:  ', bh_ID, f['PartType5']['ParticleIDs'].value, ics['PartType5']['ParticleIDs'].value
     NumPart_ThisFile[5] = 1
     ics['Header'].attrs['NumPart_ThisFile'] = NumPart_ThisFile
     for keyname in ['Masses','BH_Mass']:
       del ics['PartType5'][keyname]
       ics['PartType5'][keyname] = 10.**logMbh / (UnitMass_in_Msun/ics['Header'].attrs['HubbleParam'])
     for keyname in ['BH_Mass_AlphaDisk','BH_Mdot']:
       del ics['PartType5'][keyname]
       ics['PartType5'][keyname] = 0.

   else:

     print '  no BHs'
     NumPart_ThisFile[5] = 0
     ics['Header'].attrs['NumPart_ThisFile'] = NumPart_ThisFile
     if 'PartType5' in ics.keys():
       del ics['PartType5']


   f.close()
   ics.close()


"""
# --- create output times file:
dt_Myr = 1.0*u.Myr
nsnap = 500
cosmo = FlatLambdaCDM(H0=P['hubble']*100, Om0=P['omega_matter']) 
t_begin = cosmo.age(z_begin).to(u.Myr)
t_all = np.arange(nsnap) * dt_Myr + t_begin
z_all = np.array([])
for i in range(nsnap):
   z_all = np.append( z_all, z_at_value(cosmo.age,t_all[i]) )
a_all = 1./(1.+z_all)
# make first output very close to the ICs?
#da = a_all[1] - a_all[0]
#a_all[0] += 1e-3 * da 
filename = outdir+'/snaplist_sn%s_z%.2f_%.2fMyr.txt' % ( snapnum, z_begin, dt_Myr.value )
np.savetxt( filename, a_all, fmt='%10.8f' )
print '\n',filename
"""

print '\n\nDone!!\n'



