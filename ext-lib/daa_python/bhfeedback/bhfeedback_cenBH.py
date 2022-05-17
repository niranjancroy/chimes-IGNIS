# run as >>> execfile('bhfeedback_cenBH.py')
# python bhfeedback_cenBH.py  /simons/scratch/dangles/FIRE/bhfeedback/  h113_HR_sn214/nof_s8e1_n128  1  2   > log/cenBH.0.log 2>&1 &

import numpy as np
from scipy.stats import binned_statistic
import gadget as g
import utilities as util
import pandas as pd
import time as time
import os as os
import sys as sys
import gc as gc
import h5py as h5py
from glob import glob
from copy import deepcopy

import daa_lib as daa
from daa_constants import *

import imp
imp.reload(daa)

time_start = time.time()




#######################################
 ### DEFINE SIMULATION DIRECTORIES ###
#######################################


narg = len(sys.argv)
if narg == 5:
  basedir = sys.argv[1]
  sname = sys.argv[2]
  snap_ini = int(sys.argv[3])
  snap_end = int(sys.argv[4])
else:
  print 'syntax: python bhfeedback_cenBH.py basedir sname snap_ini snap_end'
  sys.exit()

"""
basedir = '/simons/scratch/dangles/FIRE/bhfeedback/'
sname = 'h113_HR_sn152/nof_s8e1_n128'
snap_ini = 252
snap_end = 253
"""

simdir = basedir + sname
print 'doing...', simdir


outdir = simdir + '/cenBH'
if not os.path.exists(outdir):
  try:
    os.mkdir(outdir)
  except OSError:
    try:
      os.makedirs(outdir)
    except OSError:
      print '...could not create output directory: ' + outdir




SkipIfDone = 1


ID_BAL = 4060782041
#m_wind = 3000.

#dr = 0.1
dr = 0.01
rmin = -2
rmax = 2.5 + dr
radial_bins = 10**np.arange(rmin,rmax,dr)
nbin = radial_bins.size - 1
percentiles = [ 10, 50, 90 ]
nperc = len(percentiles)
zaxis = np.array([0.,0.,1.])



vstr = {}
vstr['snapnum'] = -1
vstr['redshift'] = -1
vstr['time'] = -1

vstr['cum'] = { 'radius':radial_bins[1:],
                'Vcirc':np.zeros(nbin),
                'Ndm':np.zeros(nbin,dtype=np.int32), 'DarkMass':np.zeros(nbin),
                'Nstar':np.zeros(nbin,dtype=np.int32), 'StarMass':np.zeros(nbin), 'L_star':np.zeros([nbin,3]), 'MbulgeStar':np.zeros(nbin),
                'COMstar_pos':np.zeros([nbin,3]), 'COMstar_vel':np.zeros([nbin,3]), 'SigmaStar':np.zeros([nbin,3]),
                'Ngas':np.zeros(nbin,dtype=np.int32), 'GasMass':np.zeros(nbin), 'L_gas':np.zeros([nbin,3]), 'MbulgeGas':np.zeros(nbin),
                'COMgas_pos':np.zeros([nbin,3]), 'COMgas_vel':np.zeros([nbin,3]), 'SigmaGas':np.zeros([nbin,3]),
                'SFR':np.zeros(nbin),
                'Ngas_BAL':np.zeros(nbin,dtype=np.int32), 'GasMass_BAL':np.zeros(nbin), 'L_gas_BAL':np.zeros([nbin,3]), 'COMgas_pos_BAL':np.zeros([nbin,3]), 'COMgas_vel_BAL':np.zeros([nbin,3])
              }
 
vstr['sph'] = { 'radial_bins':radial_bins,  
                'Ndm':np.zeros(nbin,dtype=np.int32), 'DarkMass':np.zeros(nbin), 'm_dm_perc':np.zeros([nbin,nperc]),
                'Nstar':np.zeros(nbin,dtype=np.int32), 'StarMass':np.zeros(nbin), 'SigmaStar':np.zeros([nbin,3]), 'VstarPhi':np.zeros(nbin),
                'L_star':np.zeros([nbin,3]), 'm_star_perc':np.zeros([nbin,nperc]),
                'Ngas':np.zeros(nbin,dtype=np.int32), 'GasMass':np.zeros(nbin), 'SigmaGas':np.zeros([nbin,3]), 
                'L_gas':np.zeros([nbin,3]), 'm_gas_perc':np.zeros([nbin,nperc]),
                'Temp':np.zeros(nbin), 'SoundSpeed':np.zeros(nbin), 'VgasPhi':np.zeros(nbin),
                'GasMassInflow':np.zeros(nbin), 'GasInflowVel':np.zeros(nbin), 'GasInflowVel_perc':np.zeros([nbin,nperc]), 'GasInflowRate':np.zeros(nbin),
                'GasMassOutflow':np.zeros(nbin), 'GasOutflowVel':np.zeros(nbin), 'GasOutflowVel_perc':np.zeros([nbin,nperc]), 'GasOutflowRate':np.zeros(nbin),
                'SFR':np.zeros(nbin),
                'Ngas_BAL':np.zeros(nbin,dtype=np.int32), 'GasMass_BAL':np.zeros(nbin), 'm_gas_perc_BAL':np.zeros([nbin,nperc]), 'L_gas_BAL':np.zeros([nbin,3]), 'Temp_BAL':np.zeros(nbin), 'SoundSpeed_BAL':np.zeros(nbin),
                'GasMassOutflow_BAL':np.zeros(nbin), 'GasOutflowVel_BAL':np.zeros(nbin), 'GasOutflowVel_perc_BAL':np.zeros([nbin,nperc]), 'GasOutflowRate_BAL':np.zeros(nbin)
              }

cylstr = { 'radial_bins':radial_bins,
           'Nstar':np.zeros(nbin,dtype=np.int32), 'StarMass':np.zeros(nbin), 'SigmaStar':np.zeros([nbin,3]), 'L_star':np.zeros([nbin,3]), 'VstarPhi':np.zeros(nbin),
           'Ngas':np.zeros(nbin,dtype=np.int32), 'GasMass':np.zeros(nbin), 'SigmaGas':np.zeros([nbin,3]), 'L_gas':np.zeros([nbin,3]),
           'Temp':np.zeros(nbin), 'SoundSpeed':np.zeros(nbin), 'VgasPhi':np.zeros(nbin),
           'GasMassInflow':np.zeros(nbin), 'GasInflowVel':np.zeros(nbin), 'GasInflowVel_perc':np.zeros([nbin,nperc]), 'GasInflowRate':np.zeros(nbin),
           'GasMassOutflow':np.zeros(nbin), 'GasOutflowVel':np.zeros(nbin), 'GasOutflowVel_perc':np.zeros([nbin,nperc]), 'GasOutflowRate':np.zeros(nbin),
           'SFR':np.zeros(nbin)
          }
z_cyl = [ '100pc', '1kpc' ]
vstr['cyl'] = {}
for ztag in z_cyl:
   vstr['cyl'][ztag] = deepcopy(cylstr)

Mbh_test = [ 1e7, 1e8, 1e9 ]
vstr['bh'] = { 'Mass':0., 'BH_Mass':0., 'pos':np.zeros(3), 'vel':np.zeros(3),
              'RadiusOfInfluence':0., 'MassTest':np.array(Mbh_test), 'RadiusOfInfluenceTest':np.zeros(3)
             }




 ############################
  ### LOOP OVER REDSHIFT ###
 ############################

for snapnum in range(snap_ini,snap_end,1):

    outfile = outdir + '/cenBH_%04d.hdf5' % snapnum
    if SkipIfDone and os.path.exists( outfile ):
      print "\nfile exists:  %s\n" % outfile
      continue

    var = vstr.copy()
 
    var['snapnum'] = snapnum

    # --- find reference BH ---
    try:
      Pbh = g.readsnap( simdir, snapnum, 5, cosmological=1 )
    except KeyError:
      Pbh = g.readsnap( simdir, snapnum, 5, cosmological=1, skip_bh=1 )    

    if Pbh['k'] == -1:
      print '---> hey, no snapshot/BHs available: ', snapnum
      continue
    if Pbh['id'].size == 1:
      ibh = 0
    else:
      ibh = np.argmax( Pbh['mbh'] )
      print '---> hey, wrong # of BHs found in snapnum ', snapnum, '   ...', Pbh['id'].size
      #continue

    h = Pbh['hubble']
    Omega0 = Pbh['omega_matter']
    var['redshift'] = Pbh['redshift']
    var['time'] = 1e3 * util.age_of_universe( var['redshift'], h=h, Omega_M=Omega0 )  # in Myr
    hubble_factor = daa.hubble_z( var['redshift'] )

    # --- define REFERENCE FRAME ---
    cen_pos = Pbh['p'][ibh,:]
    cen_vel = Pbh['v'][ibh,:]

    var['bh']['Mass'] = Pbh['m'][ibh] * UnitMass_in_Msun
    var['bh']['BH_Mass'] = Pbh['mbh'][ibh] * UnitMass_in_Msun
    var['bh']['pos'][:] = cen_pos
    var['bh']['vel'][:] = cen_vel


    # --- DM ---
    Pdm = g.readsnap( simdir, snapnum, 1, cosmological=1 )
    r_dm = Pdm['p'][:,:] - cen_pos[:]
    R_dm = np.sqrt((r_dm*r_dm).sum(axis=1))
    indR = np.where(R_dm <= radial_bins.max())[0]
    r_dm = r_dm[indR,:]
    R_dm = R_dm[indR]
    m_dm = Pdm['m'][indR] * UnitMass_in_Msun
    

    # --- STARS ---
    Pstar = g.readsnap( simdir, snapnum, 4, cosmological=1 )
    r_star = Pstar['p'][:,:] - cen_pos[:]
    R_star = np.sqrt((r_star*r_star).sum(axis=1))
    indR = np.where(R_star <= radial_bins.max())[0]
    r_star = r_star[indR,:]
    R_star = R_star[indR]
    m_star = Pstar['m'][indR] * UnitMass_in_Msun
    v_star = Pstar['v'][indR,:] - cen_vel[:] +  hubble_factor * r_star * UnitLength_in_cm  / UnitVelocity_in_cm_per_s
    l0_star = r_star[:,1] * v_star[:,2] - r_star[:,2] * v_star[:,1]
    l1_star = r_star[:,2] * v_star[:,0] - r_star[:,0] * v_star[:,2]
    l2_star = r_star[:,0] * v_star[:,1] - r_star[:,1] * v_star[:,0]

    ind_sort = np.argsort( R_star )
    StarMassCum = np.cumsum( m_star[ind_sort] )
    ind_mass = np.where( StarMassCum >= var['bh']['Mass'] )[0]
    var['bh']['RadiusOfInfluence'] = R_star[ind_sort[ind_mass[0]]]
    for k in range(3):
       ind_mass = np.where( StarMassCum >= Mbh_test[k] )[0]
       var['bh']['RadiusOfInfluenceTest'][k] = R_star[ind_sort[ind_mass[0]]]
       ind_test = ( R_star <= var['bh']['RadiusOfInfluenceTest'][k] )
       print 'RadiusOfInfluenceTest: ', np.sum(m_star[ind_test]) / Mbh_test[k]


    # --- GAS ---
    Pgas = g.readsnap( simdir, snapnum, 0, cosmological=1 )
    r_gas = Pgas['p'][:,:] - cen_pos[:]
    R_gas = np.sqrt((r_gas*r_gas).sum(axis=1))
    indR = np.where(R_gas <= radial_bins.max())[0]
    id_gas =Pgas['id'][indR]
    r_gas = r_gas[indR,:]
    R_gas = R_gas[indR]
    m_gas = Pgas['m'][indR] * UnitMass_in_Msun
    v_gas = Pgas['v'][indR,:] - cen_vel[:] +  hubble_factor * r_gas * UnitLength_in_cm  / UnitVelocity_in_cm_per_s
    Vr_gas = (v_gas*r_gas).sum(axis=1) / R_gas
    T_gas = g.gas_temperature(Pgas['u'][indR],Pgas['ne'][indR])
    cs_gas = np.sqrt( GAMMA*GAMMA_MINUS1 * Pgas['u'][indR] )
    sfr_gas = Pgas['sfr'][indR]
    l0_gas = r_gas[:,1] * v_gas[:,2] - r_gas[:,2] * v_gas[:,1]
    l1_gas = r_gas[:,2] * v_gas[:,0] - r_gas[:,0] * v_gas[:,2]
    l2_gas = r_gas[:,0] * v_gas[:,1] - r_gas[:,1] * v_gas[:,0]



    # --- loop over all radii

    for i in range(nbin):

       print i, '/', nbin-1


       # **********************************************************************
       print '... cumulative radial properties: R < ', var['cum']['radius'][i]
       # **********************************************************************

       # --- DM
       ind_dm = np.where( R_dm < var['cum']['radius'][i] )[0]
       Ndm = ind_dm.size
       DarkMass = 0.
       if Ndm > 0:
         DarkMass = np.sum( m_dm[ind_dm] )
       var['cum']['Ndm'][i] = Ndm
       var['cum']['DarkMass'][i] = DarkMass

       # --- STARS
       ind_star = np.where( R_star < var['cum']['radius'][i] )[0]
       Nstar = ind_star.size
       Nstar_cum = Nstar
       StarMass = 0.
       MbulgeStar = 0.
       COMstar_pos = np.zeros(3)
       COMstar_vel = np.zeros(3)
       L_star = np.zeros(3)
       SigmaStar = np.zeros(3)
       if Nstar > 0:
         StarMass = np.sum( m_star[ind_star] )
         COMstar_pos = np.sum( m_star[ind_star,np.newaxis] * r_star[ind_star,:], axis=0 ) / StarMass
         COMstar_vel = np.sum( m_star[ind_star,np.newaxis] * v_star[ind_star,:], axis=0 ) / StarMass
         L_star = np.sum( np.cross(r_star[ind_star,:],v_star[ind_star,:]) * m_star[ind_star,np.newaxis], axis=0 )
         L_dir = L_star / np.linalg.norm(L_star)
         # --- rotate coordinates
         rotax = np.cross(L_dir,zaxis)
         rotangle = np.arccos(np.dot(L_dir,zaxis))
         rotmatrix = daa.rotation_matrix( rotax, rotangle )
         r_star_rot = np.tensordot( rotmatrix, r_star, axes=(1,1) ).T
         v_star_rot = np.tensordot( rotmatrix, v_star, axes=(1,1) ).T
         R_star_cyl = np.sqrt( r_star_rot[:,0]**2. + r_star_rot[:,1]**2. )
         v_star_phi = ( v_star_rot[:,0] * (-1.)*r_star_rot[:,1] + v_star_rot[:,1] * r_star_rot[:,0] ) / R_star_cyl
         # --- compute quantities in rotated frame
         vel = np.sum( m_star[ind_star,np.newaxis] * v_star_rot[ind_star,:], axis=0 ) / StarMass
         vel2 = np.sum( m_star[ind_star,np.newaxis] * v_star_rot[ind_star,:]**2, axis=0 ) / StarMass
         SigmaStar = np.sqrt( np.abs( vel2 - vel**2 ) )
         ind_bulge = np.where( v_star_phi[ind_star] < 0 )[0]
         if ind_bulge.size > 0:
           MbulgeStar = np.sum( m_star[ind_star[ind_bulge]] ) * 2.
           if ind_bulge.size > Nstar/2.:
             print '--> BULGE size =', ind_bulge.size, ' / ', Nstar 
       var['cum']['Nstar'][i] = Nstar
       var['cum']['StarMass'][i] = StarMass
       var['cum']['MbulgeStar'][i] = MbulgeStar
       var['cum']['COMstar_pos'][i,:] = COMstar_pos
       var['cum']['COMstar_vel'][i,:] = COMstar_vel
       var['cum']['L_star'][i,:] = L_star
       var['cum']['SigmaStar'][i,:] = SigmaStar

       # --- GAS
       ind_gas = np.where( (R_gas < var['cum']['radius'][i]) & (id_gas != ID_BAL) )[0]
       Ngas = ind_gas.size
       Ngas_cum = Ngas
       GasMass = 0.
       SFR = 0.
       MbulgeGas = 0.
       COMgas_pos = np.zeros(3)
       COMgas_vel = np.zeros(3)
       L_gas = np.zeros(3)
       SigmaGas = np.zeros(3)
       if Ngas > 0:
         GasMass = np.sum( m_gas[ind_gas] )
         SFR = np.sum( sfr_gas[ind_gas] )
         COMgas_pos = np.sum( m_gas[ind_gas,np.newaxis] * r_gas[ind_gas,:], axis=0 ) / GasMass
         COMgas_vel = np.sum( m_gas[ind_gas,np.newaxis] * v_gas[ind_gas,:], axis=0 ) / GasMass
         L_gas = np.sum( np.cross(r_gas[ind_gas,:],v_gas[ind_gas,:]) * m_gas[ind_gas,np.newaxis], axis=0 )
         L_dir = L_gas / np.linalg.norm(L_gas)
         # --- rotate coordinates
         rotax = np.cross(L_dir,zaxis)
         rotangle = np.arccos(np.dot(L_dir,zaxis))
         rotmatrix = daa.rotation_matrix( rotax, rotangle )
         r_gas_rot = np.tensordot( rotmatrix, r_gas, axes=(1,1) ).T
         v_gas_rot = np.tensordot( rotmatrix, v_gas, axes=(1,1) ).T
         R_gas_cyl = np.sqrt( r_gas_rot[:,0]**2. + r_gas_rot[:,1]**2. )
         v_gas_phi = ( v_gas_rot[:,0] * (-1.)*r_gas_rot[:,1] + v_gas_rot[:,1] * r_gas_rot[:,0] ) / R_gas_cyl
         Vr_gas_cyl = ( v_gas_rot[:,0]*r_gas_rot[:,0] + v_gas_rot[:,1]*r_gas_rot[:,1] ) / R_gas_cyl
         # --- compute quantities in rotated frame
         vel = np.sum( m_gas[ind_gas,np.newaxis] * v_gas_rot[ind_gas,:], axis=0 ) / GasMass
         vel2 = np.sum( m_gas[ind_gas,np.newaxis] * v_gas_rot[ind_gas,:]**2, axis=0 ) / GasMass
         SigmaGas = np.sqrt( np.abs( vel2 - vel**2 ) )
         ind_bulge = np.where( v_gas_phi[ind_gas] < 0 )[0]
         if ind_bulge.size > 0:
           MbulgeGas = np.sum( m_gas[ind_gas[ind_bulge]] ) * 2.
           if ind_bulge.size > Ngas/2.:
             print '--> BULGE size =', ind_bulge.size, ' / ', Ngas
       var['cum']['Ngas'][i] = Ngas
       var['cum']['GasMass'][i] = GasMass
       var['cum']['SFR'][i] = SFR
       var['cum']['MbulgeGas'][i] = MbulgeGas
       var['cum']['COMgas_pos'][i,:] = COMgas_pos
       var['cum']['COMgas_vel'][i,:] = COMgas_vel
       var['cum']['L_gas'][i,:] = L_gas
       var['cum']['SigmaGas'][i,:] = SigmaGas

       # --- GAS: bal wind particles
       ind_gas = np.where( (R_gas < var['cum']['radius'][i]) & (id_gas == ID_BAL) )[0]
       Ngas = ind_gas.size
       Ngas_cum_BAL = Ngas
       GasMass = 0.
       COMgas_pos = np.zeros(3)
       COMgas_vel = np.zeros(3)
       L_gas = np.zeros(3)
       if Ngas > 0:
         GasMass = np.sum( m_gas[ind_gas] )
         COMgas_pos = np.sum( m_gas[ind_gas,np.newaxis] * r_gas[ind_gas,:], axis=0 ) / GasMass
         COMgas_vel = np.sum( m_gas[ind_gas,np.newaxis] * v_gas[ind_gas,:], axis=0 ) / GasMass
         L_gas = np.sum( np.cross(r_gas[ind_gas,:],v_gas[ind_gas,:]) * m_gas[ind_gas,np.newaxis], axis=0 )
       var['cum']['Ngas_BAL'][i] = Ngas
       var['cum']['GasMass_BAL'][i] = GasMass
       var['cum']['COMgas_pos_BAL'][i,:] = COMgas_pos
       var['cum']['COMgas_vel_BAL'][i,:] = COMgas_vel
       var['cum']['L_gas_BAL'][i,:] = L_gas

       # --- ALL
       Mtot = DarkMass + StarMass + var['cum']['GasMass'][i] + var['cum']['GasMass_BAL'][i] + var['bh']['Mass']
       Vcirc = np.sqrt( G_UNIV * (Mtot*MSUN) / (var['cum']['radius'][i]*CM_PER_KPC) ) / CM_PER_KM     # km/s
       var['cum']['Vcirc'][i] = Vcirc


       # *******************************************************************************
       print '... spherical radial bins: ', radial_bins[i], ' < R < ', radial_bins[i+1]
       # *******************************************************************************

       # --- DARK MATTER
       ind_dm = np.where( (R_dm>=radial_bins[i]) & (R_dm<radial_bins[i+1]) )[0]
       Ndm = ind_dm.size
       DarkMass = 0.
       m_dm_perc = np.zeros(nperc)
       if Ndm > 0:
         DarkMass = np.sum( m_dm[ind_dm] )
         m_dm_perc[:] = np.percentile( m_dm[ind_dm], percentiles )
       var['sph']['Ndm'][i] = Ndm
       var['sph']['DarkMass'][i] = DarkMass
       var['sph']['m_dm_perc'][i,:] = m_dm_perc

       # --- STARS
       if Nstar_cum > 0:
        ind_star = np.where( (R_star>=radial_bins[i]) & (R_star<radial_bins[i+1]) )[0]
        Nstar = ind_star.size
        StarMass = 0.
        VstarPhi = 0.
        L_star = np.zeros(3)
        SigmaStar = np.zeros(3)
        m_star_perc = np.zeros(nperc)
        if Nstar > 0:
          StarMass = np.sum( m_star[ind_star] )
          L_star = np.sum( np.cross(r_star[ind_star,:],v_star[ind_star,:]) * m_star[ind_star,np.newaxis], axis=0 )
          vel = np.sum( m_star[ind_star,np.newaxis] * v_star_rot[ind_star,:], axis=0 ) / StarMass
          vel2 = np.sum( m_star[ind_star,np.newaxis] * v_star_rot[ind_star,:]**2, axis=0 ) / StarMass
          SigmaStar = np.sqrt( np.abs( vel2 - vel**2 ) )
          m_star_perc[:] = np.percentile( m_star[ind_star], percentiles )
          VstarPhi = np.sum( m_star[ind_star] * v_star_phi[ind_star] ) / StarMass
        var['sph']['Nstar'][i] = Nstar
        var['sph']['StarMass'][i] = StarMass
        var['sph']['VstarPhi'][i] = VstarPhi
        var['sph']['L_star'][i,:] = L_star
        var['sph']['SigmaStar'][i,:] = SigmaStar
        var['sph']['m_star_perc'][i,:] = m_star_perc

       # --- GAS
       if Ngas_cum > 0:
        ind_gas = np.where( (R_gas>=radial_bins[i]) & (R_gas<radial_bins[i+1]) & (id_gas != ID_BAL) )[0]
        Ngas = ind_gas.size
        GasMass = 0.
        Temp = 0.
        SoundSpeed = 0.
        VgasPhi = 0.
        SFR = 0.
        L_gas = np.zeros(3)
        SigmaGas = np.zeros(3)
        m_gas_perc = np.zeros(nperc)
        GasMassInflow = 0.
        GasInflowVel = 0.
        GasInflowVel_perc = np.zeros(nperc)
        GasInflowRate = 0.
        GasMassOutflow = 0.
        GasOutflowVel = 0.
        GasOutflowVel_perc = np.zeros(nperc)
        GasOutflowRate = 0.
        if Ngas > 0:
          GasMass = np.sum( m_gas[ind_gas] )
          Temp = np.sum( m_gas[ind_gas]*T_gas[ind_gas] ) / GasMass
          SoundSpeed = np.sum( m_gas[ind_gas]*cs_gas[ind_gas] ) / GasMass
          SFR = np.sum( sfr_gas[ind_gas] )
          L_gas = np.sum( np.cross(r_gas[ind_gas,:],v_gas[ind_gas,:]) * m_gas[ind_gas,np.newaxis], axis=0 )
          vel = np.sum( m_gas[ind_gas,np.newaxis] * v_gas_rot[ind_gas,:], axis=0 ) / GasMass
          vel2 = np.sum( m_gas[ind_gas,np.newaxis] * v_gas_rot[ind_gas,:]**2, axis=0 ) / GasMass
          SigmaGas = np.sqrt( np.abs( vel2 - vel**2 ) )
          m_gas_perc[:] = np.percentile( m_gas[ind_gas], percentiles )
          VgasPhi = np.sum( m_gas[ind_gas] * v_gas_phi[ind_gas] ) / GasMass
          dRbin = radial_bins[i+1] - radial_bins[i]
          ind_in = np.where( Vr_gas[ind_gas] < 0 )[0]
          if ind_in.size > 0:
            ind_in = ind_gas[ind_in]
            GasMassInflow = np.sum( m_gas[ind_in] )
            GasInflowVel = np.sum( Vr_gas[ind_in]*m_gas[ind_in] ) / GasMassInflow
            GasInflowVel_perc[:] = np.percentile( Vr_gas[ind_in], percentiles )
            GasInflowRate = GasMassInflow * (GasInflowVel*CM_PER_KM/CM_PER_KPC*SEC_PER_YEAR) / dRbin   # Msun/yr
          ind_out = np.where( Vr_gas[ind_gas] >= 0 )[0]
          if ind_out.size > 0:
            ind_out = ind_gas[ind_out]
            GasMassOutflow = np.sum( m_gas[ind_out] )
            GasOutflowVel = np.sum( Vr_gas[ind_out]*m_gas[ind_out] ) / GasMassOutflow
            GasOutflowVel_perc[:] = np.percentile( Vr_gas[ind_out], percentiles )
            GasOutflowRate = GasMassOutflow * (GasOutflowVel*CM_PER_KM/CM_PER_KPC*SEC_PER_YEAR) / dRbin
        var['sph']['Ngas'][i] = Ngas
        var['sph']['GasMass'][i] = GasMass
        var['sph']['Temp'][i] = Temp
        var['sph']['SoundSpeed'][i] = SoundSpeed
        var['sph']['VgasPhi'][i] = VgasPhi
        var['sph']['SFR'][i] = SFR
        var['sph']['SigmaGas'][i,:] = SigmaGas
        var['sph']['L_gas'][i,:] = L_gas
        var['sph']['m_gas_perc'][i,:] = m_gas_perc
        var['sph']['GasMassInflow'][i] = GasMassInflow
        var['sph']['GasInflowVel'][i] = GasInflowVel
        var['sph']['GasInflowVel_perc'][i,:] = GasInflowVel_perc
        var['sph']['GasInflowRate'][i] = GasInflowRate
        var['sph']['GasMassOutflow'][i] = GasMassOutflow
        var['sph']['GasOutflowVel'][i] = GasOutflowVel
        var['sph']['GasOutflowVel_perc'][i,:] = GasOutflowVel_perc
        var['sph']['GasOutflowRate'][i] = GasOutflowRate

       # --- GAS: bal wind particles
       if Ngas_cum_BAL > 0:
        ind_gas = np.where( (R_gas>=radial_bins[i]) & (R_gas<radial_bins[i+1]) & (id_gas == ID_BAL) )[0]
        Ngas = ind_gas.size
        GasMass = 0.
        Temp = 0.
        SoundSpeed = 0.
        L_gas = np.zeros(3)
        m_gas_perc = np.zeros(nperc)
        GasMassOutflow = 0.
        GasOutflowVel = 0.
        GasOutflowVel_perc = np.zeros(nperc)
        GasOutflowRate = 0.
        if Ngas > 0:
          GasMass = np.sum( m_gas[ind_gas] )
          Temp = np.sum( m_gas[ind_gas]*T_gas[ind_gas] ) / GasMass
          SoundSpeed = np.sum( m_gas[ind_gas]*cs_gas[ind_gas] ) / GasMass
          L_gas = np.sum( np.cross(r_gas[ind_gas,:],v_gas[ind_gas,:]) * m_gas[ind_gas,np.newaxis], axis=0 )
          m_gas_perc[:] = np.percentile( m_gas[ind_gas], percentiles )
          dRbin = radial_bins[i+1] - radial_bins[i]
          ind_out = np.where( Vr_gas[ind_gas] >= 0 )[0]
          if ind_out.size > 0:
            ind_out = ind_gas[ind_out]
            GasMassOutflow = np.sum( m_gas[ind_out] )
            GasOutflowVel = np.sum( Vr_gas[ind_out]*m_gas[ind_out] ) / GasMassOutflow
            GasOutflowVel_perc[:] = np.percentile( Vr_gas[ind_out], percentiles )
            GasOutflowRate = GasMassOutflow * (GasOutflowVel*CM_PER_KM/CM_PER_KPC*SEC_PER_YEAR) / dRbin
        var['sph']['Ngas_BAL'][i] = Ngas
        var['sph']['GasMass_BAL'][i] = GasMass
        var['sph']['Temp_BAL'][i] = Temp
        var['sph']['SoundSpeed_BAL'][i] = SoundSpeed
        var['sph']['L_gas_BAL'][i,:] = L_gas
        var['sph']['m_gas_perc_BAL'][i,:] = m_gas_perc
        var['sph']['GasMassOutflow_BAL'][i] = GasMassOutflow
        var['sph']['GasOutflowVel_BAL'][i] = GasOutflowVel
        var['sph']['GasOutflowVel_perc_BAL'][i,:] = GasOutflowVel_perc
        var['sph']['GasOutflowRate_BAL'][i] = GasOutflowRate        


       for ztag, z_max in zip( ['100pc', '1kpc'], [0.1, 1.] ):

          # **********************************************************************************************************
          print '... cylindrical radial bins: ', radial_bins[i], ' < R_cyl < ', radial_bins[i+1], ',   |z| < ', z_max
          # **********************************************************************************************************

          # --- STARS
          if Nstar_cum > 0:
           ind_star = np.where( (R_star_cyl>=radial_bins[i]) & (R_star_cyl<radial_bins[i+1]) & (np.abs(r_star_rot[:,2])<z_max) )[0]
           Nstar = ind_star.size
           StarMass = 0.
           VstarPhi = 0.
           L_star = np.zeros(3)
           SigmaStar = np.zeros(3)
           if Nstar > 0:
             StarMass = np.sum( m_star[ind_star] )
             L_star = np.sum( np.cross(r_star[ind_star,:],v_star[ind_star,:]) * m_star[ind_star,np.newaxis], axis=0 )
             vel = np.sum( m_star[ind_star,np.newaxis] * v_star_rot[ind_star,:], axis=0 ) / StarMass
             vel2 = np.sum( m_star[ind_star,np.newaxis] * v_star_rot[ind_star,:]**2, axis=0 ) / StarMass
             SigmaStar = np.sqrt( np.abs( vel2 - vel**2 ) )
             VstarPhi = np.sum( m_star[ind_star] * v_star_phi[ind_star] ) / StarMass
           var['cyl'][ztag]['Nstar'][i] = Nstar
           var['cyl'][ztag]['StarMass'][i] = StarMass
           var['cyl'][ztag]['VstarPhi'][i] = VstarPhi
           var['cyl'][ztag]['L_star'][i,:] = L_star
           var['cyl'][ztag]['SigmaStar'][i,:] = SigmaStar

          # --- GAS
          if Ngas_cum > 0:
           ind_gas = np.where( (R_gas_cyl>=radial_bins[i]) & (R_gas_cyl<radial_bins[i+1]) & (np.abs(r_gas_rot[:,2])<z_max) & (id_gas != ID_BAL) )[0]
           Ngas = ind_gas.size
           GasMass = 0.
           Temp = 0.
           SoundSpeed = 0.
           VgasPhi = 0.
           SFR = 0.
           L_gas = np.zeros(3)
           SigmaGas = np.zeros(3)
           GasMassInflow = 0.
           GasInflowVel = 0.
           GasInflowVel_perc = np.zeros(nperc)
           GasInflowRate = 0.
           GasMassOutflow = 0.
           GasOutflowVel = 0.
           GasOutflowVel_perc = np.zeros(nperc)
           GasOutflowRate = 0.
           if Ngas > 0:
             GasMass = np.sum( m_gas[ind_gas] )
             Temp = np.sum( m_gas[ind_gas]*T_gas[ind_gas] ) / GasMass
             SoundSpeed = np.sum( m_gas[ind_gas]*cs_gas[ind_gas] ) / GasMass
             SFR = np.sum( sfr_gas[ind_gas] )
             L_gas = np.sum( np.cross(r_gas[ind_gas,:],v_gas[ind_gas,:]) * m_gas[ind_gas,np.newaxis], axis=0 )
             vel = np.sum( m_gas[ind_gas,np.newaxis] * v_gas_rot[ind_gas,:], axis=0 ) / GasMass
             vel2 = np.sum( m_gas[ind_gas,np.newaxis] * v_gas_rot[ind_gas,:]**2, axis=0 ) / GasMass
             SigmaGas = np.sqrt( np.abs( vel2 - vel**2 ) )
             VgasPhi = np.sum( m_gas[ind_gas] * v_gas_phi[ind_gas] ) / GasMass
             dRbin = radial_bins[i+1] - radial_bins[i]
             ind_in = np.where( Vr_gas_cyl[ind_gas] < 0 )[0]
             if ind_in.size > 0:
               ind_in = ind_gas[ind_in]
               GasMassInflow = np.sum( m_gas[ind_in] )
               GasInflowVel = np.sum( Vr_gas_cyl[ind_in]*m_gas[ind_in] ) / GasMassInflow
               GasInflowVel_perc[:] = np.percentile( Vr_gas_cyl[ind_in], percentiles )
               GasInflowRate = GasMassInflow * (GasInflowVel*CM_PER_KM/CM_PER_KPC*SEC_PER_YEAR) / dRbin   # Msun/yr
             ind_out = np.where( Vr_gas_cyl[ind_gas] >= 0 )[0]
             if ind_out.size > 0:
               ind_out = ind_gas[ind_out]
               GasMassOutflow = np.sum( m_gas[ind_out] )
               GasOutflowVel = np.sum( Vr_gas_cyl[ind_out]*m_gas[ind_out] ) / GasMassOutflow
               GasOutflowVel_perc[:] = np.percentile( Vr_gas_cyl[ind_out], percentiles )
               GasOutflowRate = GasMassOutflow * (GasOutflowVel*CM_PER_KM/CM_PER_KPC*SEC_PER_YEAR) / dRbin
           var['cyl'][ztag]['Ngas'][i] = Ngas
           var['cyl'][ztag]['GasMass'][i] = GasMass
           var['cyl'][ztag]['Temp'][i] = Temp
           var['cyl'][ztag]['SoundSpeed'][i] = SoundSpeed
           var['cyl'][ztag]['VgasPhi'][i] = VgasPhi
           var['cyl'][ztag]['SFR'][i] = SFR
           var['cyl'][ztag]['SigmaGas'][i,:] = SigmaGas
           var['cyl'][ztag]['L_gas'][i,:] = L_gas
           var['cyl'][ztag]['GasMassInflow'][i] = GasMassInflow
           var['cyl'][ztag]['GasInflowVel'][i] = GasInflowVel
           var['cyl'][ztag]['GasInflowVel_perc'][i,:] = GasInflowVel_perc
           var['cyl'][ztag]['GasInflowRate'][i] = GasInflowRate
           var['cyl'][ztag]['GasMassOutflow'][i] = GasMassOutflow
           var['cyl'][ztag]['GasOutflowVel'][i] = GasOutflowVel
           var['cyl'][ztag]['GasOutflowVel_perc'][i,:] = GasOutflowVel_perc
           var['cyl'][ztag]['GasOutflowRate'][i] = GasOutflowRate  


    print snapnum, ':  redshift =', var['redshift']
    sys.stdout.flush()


    ############################
     ### WRITE DATA TO FILE ###
    ############################

    print '\n... writing', outfile
    if os.path.isfile(outfile ):
      os.remove(outfile)
    file = h5py.File(outfile, 'w')
    for grname in var.keys():
       if type(var[grname]) is dict:
         group = file.create_group( grname )
         for keyname in var[grname].keys():
            if type(var[grname][keyname]) is dict:
              group2 = group.create_group( keyname )
              for keyname2 in var[grname][keyname].keys():
                 group2.create_dataset(keyname2, data=var[grname][keyname][keyname2])
            else:
              group.create_dataset(keyname, data=var[grname][keyname])
       else:
         file.create_dataset(grname, data=var[grname])
    file.close()




print '\n\nDone!!'

