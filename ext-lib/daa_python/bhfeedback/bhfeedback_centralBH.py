# run as >>> execfile('bhfeedback_centralBH.py')
# python bhfeedback_centralBH.py > log/cenBH.0.log 2>&1 &
# python bhfeedback_centralBH.py bal_s9e1v30f05_n256 > log/cenBH.0.log 2>&1 &

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
if narg == 2:
  sname = sys.argv[1]
  sim_list = [sname]
else:
  print 'syntax: python bhfeedback_centralBH.py sname'
  sys.exit()


"""
sim_list = [ 'spawn_s9e1v30f05_w1m1000_new' ]
#isim = 15
#sim_list = [ sim_list[isim] ]
"""

nsim = len(sim_list)

basedir = '/simons/scratch/dangles/FIRE/bhfeedback/'
#ICdir = 'h113_HR_sn214/'

for sname in sim_list:

 #simdir = basedir + ICdir + sname
 simdir = basedir + sname
 print 'doing...', simdir

 snapdir = glob( simdir + '/snapdir*')
 nsnap = len(snapdir)
 snap_list = np.zeros(nsnap,dtype='int')
 for i in range(nsnap):
    snap_list[i] = int( snapdir[i][ snapdir[i].find('_', -5, -1)+1: ] )
 snap_list.sort()
 print snap_list

 var = {}
 var['snapnum'] = np.zeros(nsnap,dtype=np.int32)
 var['redshift'] = np.zeros(nsnap)
 var['time'] = np.zeros(nsnap)

 dr = 0.1
 rmin = -2
 rmax = 2 + dr
 radial_bins = 10**np.arange(rmin,rmax,dr)
 nbin = radial_bins.size - 1
 percentiles = [ 10, 50, 90 ]
 nperc = len(percentiles)
 zaxis = np.array([0.,0.,1.])

 var['cum'] = { 'radius':radial_bins[1:],
                'Vcirc':np.zeros([nsnap,nbin]),
                'Ndm':np.zeros([nsnap,nbin],dtype=np.int32), 'DarkMass':np.zeros([nsnap,nbin]),
                'Nstar':np.zeros([nsnap,nbin],dtype=np.int32), 'StarMass':np.zeros([nsnap,nbin]), 'L_star':np.zeros([nsnap,nbin,3]), 'MbulgeStar':np.zeros([nsnap,nbin]),
                'COMstar_pos':np.zeros([nsnap,nbin,3]), 'COMstar_vel':np.zeros([nsnap,nbin,3]), 'SigmaStar':np.zeros([nsnap,nbin,3]),
                'Ngas':np.zeros([nsnap,nbin],dtype=np.int32), 'GasMass':np.zeros([nsnap,nbin]), 'L_gas':np.zeros([nsnap,nbin,3]), 'MbulgeGas':np.zeros([nsnap,nbin]),
                'COMgas_pos':np.zeros([nsnap,nbin,3]), 'COMgas_vel':np.zeros([nsnap,nbin,3]), 'SigmaGas':np.zeros([nsnap,nbin,3]),
                'SFR':np.zeros([nsnap,nbin])
              }
 
 var['sph'] = { 'radial_bins':radial_bins,  
                'Ndm':np.zeros([nsnap,nbin],dtype=np.int32), 'DarkMass':np.zeros([nsnap,nbin]), 'm_dm_perc':np.zeros([nsnap,nbin,nperc]),
                'Nstar':np.zeros([nsnap,nbin],dtype=np.int32), 'StarMass':np.zeros([nsnap,nbin]), 'SigmaStar':np.zeros([nsnap,nbin,3]), 'VstarPhi':np.zeros([nsnap,nbin]),
                'L_star':np.zeros([nsnap,nbin,3]), 'm_star_perc':np.zeros([nsnap,nbin,nperc]),
                'Ngas':np.zeros([nsnap,nbin],dtype=np.int32), 'GasMass':np.zeros([nsnap,nbin]), 'SigmaGas':np.zeros([nsnap,nbin,3]), 
                'L_gas':np.zeros([nsnap,nbin,3]), 'm_gas_perc':np.zeros([nsnap,nbin,nperc]),
                'Temp':np.zeros([nsnap,nbin]), 'SoundSpeed':np.zeros([nsnap,nbin]), 'VgasPhi':np.zeros([nsnap,nbin]),
                'GasMassInflow':np.zeros([nsnap,nbin]), 'GasInflowVel':np.zeros([nsnap,nbin]), 'GasInflowVel_perc':np.zeros([nsnap,nbin,nperc]), 'GasInflowRate':np.zeros([nsnap,nbin]),
                'GasMassOutflow':np.zeros([nsnap,nbin]), 'GasOutflowVel':np.zeros([nsnap,nbin]), 'GasOutflowVel_perc':np.zeros([nsnap,nbin,nperc]), 'GasOutflowRate':np.zeros([nsnap,nbin]),
                'SFR':np.zeros([nsnap,nbin])
              }

 cylstr = { 'radial_bins':radial_bins,
            'Nstar':np.zeros([nsnap,nbin],dtype=np.int32), 'StarMass':np.zeros([nsnap,nbin]), 'SigmaStar':np.zeros([nsnap,nbin,3]), 'L_star':np.zeros([nsnap,nbin,3]), 'VstarPhi':np.zeros([nsnap,nbin]),
            'Ngas':np.zeros([nsnap,nbin],dtype=np.int32), 'GasMass':np.zeros([nsnap,nbin]), 'SigmaGas':np.zeros([nsnap,nbin,3]), 'L_gas':np.zeros([nsnap,nbin,3]),
            'Temp':np.zeros([nsnap,nbin]), 'SoundSpeed':np.zeros([nsnap,nbin]), 'VgasPhi':np.zeros([nsnap,nbin]),
            'GasMassInflow':np.zeros([nsnap,nbin]), 'GasInflowVel':np.zeros([nsnap,nbin]), 'GasInflowVel_perc':np.zeros([nsnap,nbin,nperc]), 'GasInflowRate':np.zeros([nsnap,nbin]),
            'GasMassOutflow':np.zeros([nsnap,nbin]), 'GasOutflowVel':np.zeros([nsnap,nbin]), 'GasOutflowVel_perc':np.zeros([nsnap,nbin,nperc]), 'GasOutflowRate':np.zeros([nsnap,nbin]),
            'SFR':np.zeros([nsnap,nbin])
          }
 z_cyl = [ '10pc', '100pc', '1kpc', '3th' ]
 var['cyl'] = {}
 for ztag in z_cyl:
    var['cyl'][ztag] = deepcopy(cylstr)

 Mbh_test = [ 1e7, 1e8, 1e9 ]
 var['bh'] = { 'Mass':np.zeros(nsnap), 'BH_Mass':np.zeros(nsnap), 'pos':np.zeros([nsnap,3]), 'vel':np.zeros([nsnap,3]),
               'RadiusOfInfluence':np.zeros(nsnap), 'MassTest':np.array(Mbh_test), 'RadiusOfInfluenceTest':np.zeros([nsnap,3])
             }



 ############################
  ### LOOP OVER REDSHIFT ###
 ############################

 for n in range(nsnap):
 #for n in [0,nsnap-1]:
 #for n in [0]:

    snapnum = snap_list[n]
    var['snapnum'][n] = snapnum

    # --- find reference BH ---
    try:
      Pbh = g.readsnap( simdir, snapnum, 5, cosmological=1 )
    except KeyError:
      Pbh = g.readsnap( simdir, snapnum, 5, cosmological=1, skip_bh=1 )    

    if Pbh['k'] == -1:
      print '---> hey, no snapshot/BHs available: ', snapnum
      continue
    if Pbh['id'].size != 1:
      print '---> hey, wrong # of BHs found in snapnum ', snapnum, '   ...', Pbh['id'].size
      continue

    h = Pbh['hubble']
    Omega0 = Pbh['omega_matter']
    var['redshift'][n] = Pbh['redshift']
    var['time'][n] = 1e3 * util.age_of_universe( var['redshift'][n], h=h, Omega_M=Omega0 )  # in Myr
    hubble_factor = daa.hubble_z( var['redshift'][n] )

    # --- define REFERENCE FRAME ---
    cen_pos = Pbh['p'][0,:]
    cen_vel = Pbh['v'][0,:]

    var['bh']['Mass'][n] = Pbh['m'][0] * UnitMass_in_Msun
    var['bh']['BH_Mass'][n] = Pbh['mbh'][0] * UnitMass_in_Msun
    var['bh']['pos'][n,:] = cen_pos
    var['bh']['vel'][n,:] = cen_vel


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
    ind_mass = np.where( StarMassCum >= var['bh']['Mass'][n] )[0]
    var['bh']['RadiusOfInfluence'][n] = R_star[ind_sort[ind_mass[0]]]
    for k in range(3):
       ind_mass = np.where( StarMassCum >= Mbh_test[k] )[0]
       var['bh']['RadiusOfInfluenceTest'][n,k] = R_star[ind_sort[ind_mass[0]]]
       ind_test = ( R_star <= var['bh']['RadiusOfInfluenceTest'][n,k] )
       print 'RadiusOfInfluenceTest: ', np.sum(m_star[ind_test]) / Mbh_test[k]


    # --- GAS ---
    Pgas = g.readsnap( simdir, snapnum, 0, cosmological=1 )
    r_gas = Pgas['p'][:,:] - cen_pos[:]
    R_gas = np.sqrt((r_gas*r_gas).sum(axis=1))
    indR = np.where(R_gas <= radial_bins.max())[0]
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
       var['cum']['Ndm'][n,i] = Ndm
       var['cum']['DarkMass'][n,i] = DarkMass

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
       var['cum']['Nstar'][n,i] = Nstar
       var['cum']['StarMass'][n,i] = StarMass
       var['cum']['MbulgeStar'][n,i] = MbulgeStar
       var['cum']['COMstar_pos'][n,i,:] = COMstar_pos
       var['cum']['COMstar_vel'][n,i,:] = COMstar_vel
       var['cum']['L_star'][n,i,:] = L_star
       var['cum']['SigmaStar'][n,i,:] = SigmaStar

       # --- GAS
       ind_gas = np.where( R_gas < var['cum']['radius'][i] )[0]
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
       var['cum']['Ngas'][n,i] = Ngas
       var['cum']['GasMass'][n,i] = GasMass
       var['cum']['SFR'][n,i] = SFR
       var['cum']['MbulgeGas'][n,i] = MbulgeGas
       var['cum']['COMgas_pos'][n,i,:] = COMgas_pos
       var['cum']['COMgas_vel'][n,i,:] = COMgas_vel
       var['cum']['L_gas'][n,i,:] = L_gas
       var['cum']['SigmaGas'][n,i,:] = SigmaGas

       # --- ALL
       Mtot = DarkMass + StarMass + GasMass + var['bh']['Mass'][n]
       Vcirc = np.sqrt( G_UNIV * (Mtot*MSUN) / (var['cum']['radius'][i]*CM_PER_KPC) ) / CM_PER_KM     # km/s
       var['cum']['Vcirc'][n,i] = Vcirc


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
       var['sph']['Ndm'][n,i] = Ndm
       var['sph']['DarkMass'][n,i] = DarkMass
       var['sph']['m_dm_perc'][n,i,:] = m_dm_perc

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
          """
          # --- just testing...
          L_star_tmp = np.array( [ np.sum(l0_star[ind_star]*m_star[ind_star]),
                               np.sum(l1_star[ind_star]*m_star[ind_star]),
                               np.sum(l2_star[ind_star]*m_star[ind_star]) ] )
          print 'L_star_tmp/L_star = ', L_star_tmp/L_star
          """
        var['sph']['Nstar'][n,i] = Nstar
        var['sph']['StarMass'][n,i] = StarMass
        var['sph']['VstarPhi'][n,i] = VstarPhi
        var['sph']['L_star'][n,i,:] = L_star
        var['sph']['SigmaStar'][n,i,:] = SigmaStar
        var['sph']['m_star_perc'][n,i,:] = m_star_perc


       # --- GAS
       if Ngas_cum > 0:
        ind_gas = np.where( (R_gas>=radial_bins[i]) & (R_gas<radial_bins[i+1]) )[0]
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
          """
          # --- just testing...
          L_gas_tmp = var['cum']['L_gas'][n,i,:]
          L_dir_tmp = L_gas_tmp / np.linalg.norm(L_gas_tmp)
          v_gas_z = v_gas[ind_gas,0] * L_dir_tmp[0] + v_gas[ind_gas,1] * L_dir_tmp[1] + v_gas[ind_gas,2] * L_dir_tmp[2]
          vel = np.sum( m_gas[ind_gas] * v_gas_z ) / GasMass
          vel2 = np.sum( m_gas[ind_gas] * v_gas_z**2 ) / GasMass
          SigmaGasZ = np.sqrt( np.abs( vel2 - vel**2 ) )
          print 'SigmaGasZ/SigmaGas[2] = ', SigmaGasZ/SigmaGas[2]
          phi_dir = np.cross( L_dir_tmp, r_gas, axisb=1 )
          phi_dir /= np.linalg.norm(phi_dir,axis=1)[:,np.newaxis]
          v_gas_phi_tmp = v_gas[ind_gas,0] * phi_dir[ind_gas,0] + v_gas[ind_gas,1] * phi_dir[ind_gas,1] + v_gas[ind_gas,2] * phi_dir[ind_gas,2]
          VgasPhi_tmp = np.sum( m_gas[ind_gas]*v_gas_phi_tmp ) / GasMass
          print 'VgasPhi_tmp/VgasPhi = ', VgasPhi_tmp/VgasPhi
          Vr_gas_rot_tmp = (v_gas_rot*r_gas_rot).sum(axis=1) / R_gas
          ind_out = np.where( Vr_gas_rot_tmp[ind_gas] >= 0 )[0]
          if ind_out.size > 0:
            ind_out = ind_gas[ind_out]
            GasMassOutflow_tmp = np.sum( m_gas[ind_out] )
            GasOutflowVel_tmp = np.sum( Vr_gas_rot_tmp[ind_out]*m_gas[ind_out] ) / GasMassOutflow_tmp
            print 'GasMassOutflow_tmp/GasMassOutflow =', GasMassOutflow_tmp/GasMassOutflow, '  GasOutflowVel_tmp/GasOutflowVel =', GasOutflowVel_tmp/GasOutflowVel
          """
        var['sph']['Ngas'][n,i] = Ngas
        var['sph']['GasMass'][n,i] = GasMass
        var['sph']['Temp'][n,i] = Temp
        var['sph']['SoundSpeed'][n,i] = SoundSpeed
        var['sph']['VgasPhi'][n,i] = VgasPhi
        var['sph']['SFR'][n,i] = SFR
        var['sph']['SigmaGas'][n,i,:] = SigmaGas
        var['sph']['L_gas'][n,i,:] = L_gas
        var['sph']['m_gas_perc'][n,i,:] = m_gas_perc
        var['sph']['GasMassInflow'][n,i] = GasMassInflow
        var['sph']['GasInflowVel'][n,i] = GasInflowVel
        var['sph']['GasInflowVel_perc'][n,i,:] = GasInflowVel_perc
        var['sph']['GasInflowRate'][n,i] = GasInflowRate
        var['sph']['GasMassOutflow'][n,i] = GasMassOutflow
        var['sph']['GasOutflowVel'][n,i] = GasOutflowVel
        var['sph']['GasOutflowVel_perc'][n,i,:] = GasOutflowVel_perc
        var['sph']['GasOutflowRate'][n,i] = GasOutflowRate



       for ztag, z_max in zip( ['10pc', '100pc', '1kpc', '3th'], [0.01, 0.1, 1., radial_bins[i+1]/3.] ):

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
           var['cyl'][ztag]['Nstar'][n,i] = Nstar
           var['cyl'][ztag]['StarMass'][n,i] = StarMass
           var['cyl'][ztag]['VstarPhi'][n,i] = VstarPhi
           var['cyl'][ztag]['L_star'][n,i,:] = L_star
           var['cyl'][ztag]['SigmaStar'][n,i,:] = SigmaStar

          # --- GAS
          if Ngas_cum > 0:
           ind_gas = np.where( (R_gas_cyl>=radial_bins[i]) & (R_gas_cyl<radial_bins[i+1]) & (np.abs(r_gas_rot[:,2])<z_max) )[0]
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
           var['cyl'][ztag]['Ngas'][n,i] = Ngas
           var['cyl'][ztag]['GasMass'][n,i] = GasMass
           var['cyl'][ztag]['Temp'][n,i] = Temp
           var['cyl'][ztag]['SoundSpeed'][n,i] = SoundSpeed
           var['cyl'][ztag]['VgasPhi'][n,i] = VgasPhi
           var['cyl'][ztag]['SFR'][n,i] = SFR
           var['cyl'][ztag]['SigmaGas'][n,i,:] = SigmaGas
           var['cyl'][ztag]['L_gas'][n,i,:] = L_gas
           var['cyl'][ztag]['GasMassInflow'][n,i] = GasMassInflow
           var['cyl'][ztag]['GasInflowVel'][n,i] = GasInflowVel
           var['cyl'][ztag]['GasInflowVel_perc'][n,i,:] = GasInflowVel_perc
           var['cyl'][ztag]['GasInflowRate'][n,i] = GasInflowRate
           var['cyl'][ztag]['GasMassOutflow'][n,i] = GasMassOutflow
           var['cyl'][ztag]['GasOutflowVel'][n,i] = GasOutflowVel
           var['cyl'][ztag]['GasOutflowVel_perc'][n,i,:] = GasOutflowVel_perc
           var['cyl'][ztag]['GasOutflowRate'][n,i] = GasOutflowRate  



    print snapnum, ':  redshift =', var['redshift'][n]
    sys.stdout.flush()


 ############################
  ### WRITE DATA TO FILE ###
 ############################

 #outfile = simdir + '/bhdata.hdf5'
 #outfile = simdir + '/bh_rbins.hdf5'
 outfile = simdir + '/centralBH.hdf5'
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

