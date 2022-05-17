# run as >>> execfile('bhfeedback_centralBH.py')
# python bhfeedback_centralBH.py > log/cenBH.0.log 2>&1 &

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

basedir = '/simons/scratch/dangles/FIRE/bhfeedback/'
ICdir = 'h113_HR_sn214/'
sim_list = [ 'nof_s8e1_n256',
             'bal_s9e1v30f05_n256',
             'bal_s8e1v30f05_n256',
             'bal_s7e1v30f05_n256',
             'kick_s8e1v15000p2_n256',
             'kick_s8e1v1500p20_n256',
             'kick_s7e1v15000p2_n256',
             'kick_s7e1v1500p20_n256',
             'spawn_s9e1v30f05_w1m100',
             'spawn_s8e1v30f05_w1m100',
             'spawn_s7e1v30f05_w1m100',
             'spawn_s9e1v30f05_w1m10',
             'spawn_s8e1v30f05_w1m10',
             'spawn_s7e1v30f05_w1m10' ]


isim = 13
sim_list = [ sim_list[isim] ]

nsim = len(sim_list)


for sname in sim_list:

 simdir = basedir + ICdir + sname
 print 'doing...', simdir

 snapdir = glob( simdir + '/snapdir*')
 snapdir.sort()
 nsnap = len(snapdir)


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

 var['prf'] = { 'radial_bins':radial_bins,  
                'Ndm':np.zeros([nsnap,nbin],dtype=np.int32), 'DarkMass':np.zeros([nsnap,nbin]),
                'Nstar':np.zeros([nsnap,nbin],dtype=np.int32), 'StarMass':np.zeros([nsnap,nbin]), 'MbulgeStar':np.zeros([nsnap,nbin]), 'SigmaStar':np.zeros([nsnap,nbin,3]), 
                'L_star':np.zeros([nsnap,nbin,3]), 'm_star_perc':np.zeros([nsnap,nbin,nperc]),
                'Ngas':np.zeros([nsnap,nbin],dtype=np.int32), 'GasMass':np.zeros([nsnap,nbin]), 'MbulgeGas':np.zeros([nsnap,nbin]), 'SigmaGas':np.zeros([nsnap,nbin,3]), 'L_gas':np.zeros([nsnap,nbin,3]),
                'Temp':np.zeros([nsnap,nbin]), 'SoundSpeed':np.zeros([nsnap,nbin]), 'SigmaGasZ':np.zeros([nsnap,nbin]), 'VgasPhi':np.zeros([nsnap,nbin]), 'm_gas_perc':np.zeros([nsnap,nbin,nperc]), 
                'GasMassInflow':np.zeros([nsnap,nbin]), 'GasInflowVel':np.zeros([nsnap,nbin]), 'GasInflowVel_perc':np.zeros([nsnap,nbin,nperc]),
                'GasMassOutflow':np.zeros([nsnap,nbin]), 'GasOutflowVel':np.zeros([nsnap,nbin]), 'GasOutflowVel_perc':np.zeros([nsnap,nbin,nperc]),
                'SFR':np.zeros([nsnap,nbin])
              }


 ############################
  ### LOOP OVER REDSHIFT ###
 ############################

 for n in range(nsnap):
 #for n in [0,nsnap-1]:

    snapnum = int(snapdir[n][-3:])
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


    # --- DM ---
    Pdm = g.readsnap( simdir, snapnum, 1, cosmological=1 )
    r_dm = Pdm['p'][:,:] - cen_pos[:]
    R_dm = np.sqrt((r_dm*r_dm).sum(axis=1))
    m_dm = Pdm['m'][:] * UnitMass_in_Msun
    

    # --- STARS ---
    Pstar = g.readsnap( simdir, snapnum, 4, cosmological=1 )
    r_star = Pstar['p'][:,:] - cen_pos[:]
    R_star = np.sqrt((r_star*r_star).sum(axis=1))
    m_star = Pstar['m'][:] * UnitMass_in_Msun
    v_star = Pstar['v'][:,:] - cen_vel[:] +  hubble_factor * r_star * UnitLength_in_cm  / UnitVelocity_in_cm_per_s
    l0_star = r_star[:,1] * v_star[:,2] - r_star[:,2] * v_star[:,1]
    l1_star = r_star[:,2] * v_star[:,0] - r_star[:,0] * v_star[:,2]
    l2_star = r_star[:,0] * v_star[:,1] - r_star[:,1] * v_star[:,0]
 

    # --- GAS ---
    Pgas = g.readsnap( simdir, snapnum, 0, cosmological=1 )
    m_gas = Pgas['m'][:] * UnitMass_in_Msun
    r_gas = Pgas['p'][:,:] - cen_pos[:]
    R_gas = np.sqrt((r_gas*r_gas).sum(axis=1))
    v_gas = Pgas['v'][:,:] - cen_vel[:] +  hubble_factor * r_gas * UnitLength_in_cm  / UnitVelocity_in_cm_per_s
    Vr_gas = (v_gas*r_gas).sum(axis=1) / R_gas
    T_gas = g.gas_temperature(Pgas['u'],Pgas['ne'])
    cs_gas = np.sqrt( GAMMA*GAMMA_MINUS1 * Pgas['u'] )
    sfr_gas = Pgas['sfr'][:]
    l0_gas = r_gas[:,1] * v_gas[:,2] - r_gas[:,2] * v_gas[:,1]
    l1_gas = r_gas[:,2] * v_gas[:,0] - r_gas[:,0] * v_gas[:,2]
    l2_gas = r_gas[:,0] * v_gas[:,1] - r_gas[:,1] * v_gas[:,0]


    # --- compute properties within each radial bin

    for i in range(nbin):

       print i, '/', nbin-1

       # --- DARK MATTER
       ind_dm = np.where( (R_dm>=radial_bins[i]) & (R_dm<radial_bins[i+1]) )[0]

       Ndm = ind_dm.size
       DarkMass = 0.
       if Ndm > 0:
         DarkMass = np.sum( m_dm[ind_dm] )

       # --- STARS
       ind_star = np.where( (R_star>=radial_bins[i]) & (R_star<radial_bins[i+1]) )[0]
    
       Nstar = ind_star.size
       StarMass = 0.
       MbulgeStar = 0.
       L_star = np.zeros(3)
       SigmaStar = np.zeros(3)
       m_star_perc = np.zeros(nperc)
       if Nstar > 0:
         StarMass = np.sum( m_star[ind_star] )
         L_star = np.array( [ np.sum(l0_star[ind_star]*m_star[ind_star]),
                              np.sum(l1_star[ind_star]*m_star[ind_star]),
                              np.sum(l2_star[ind_star]*m_star[ind_star]) ] )
         L_dir = L_star / np.sqrt( np.dot(L_star,L_star) )
         j_star = l0_star[ind_star] * L_dir[0] + l1_star[ind_star] * L_dir[1] + l2_star[ind_star] * L_dir[2]
         ind_bulge = np.where( j_star < 0 )[0]
         MbulgeStar = np.sum( m_star[ind_star[ind_bulge]] ) * 2.
         for k in range(3):
            vel = np.sum( m_star[ind_star] * v_star[ind_star,k] ) / StarMass
            vel2 = np.sum( m_star[ind_star] * v_star[ind_star,k]**2 ) / StarMass
            SigmaStar[k] = np.sqrt( np.abs( vel2 - vel**2 ) )
         m_star_perc[:] = np.percentile( m_star[ind_star], percentiles )

       # --- GAS
       ind_gas = np.where( (R_gas>=radial_bins[i]) & (R_gas<radial_bins[i+1]) )[0]

       Ngas = ind_gas.size
       GasMass = 0.
       MbulgeGas = 0.
       Temp = 0.
       SoundSpeed = 0.
       SigmaGasZ = 0.
       VgasPhi = 0.
       SFR = 0.
       L_gas = np.zeros(3)
       SigmaGas = np.zeros(3)
       m_gas_perc = np.zeros(nperc)
       GasMassInflow = 0.
       GasInflowVel = 0.
       GasInflowVel_perc = np.zeros(nperc)
       GasMassOutflow = 0.
       GasOutflowVel = 0.
       GasOutflowVel_perc = np.zeros(nperc)
       if Ngas > 0:
         GasMass = np.sum( m_gas[ind_gas] )
         Temp = np.sum( m_gas[ind_gas]*T_gas[ind_gas] ) / GasMass
         SoundSpeed = np.sum( m_gas[ind_gas]*cs_gas[ind_gas] ) / GasMass
         SFR = np.sum( sfr_gas[ind_gas] )
         L_gas = np.array( [ np.sum(l0_gas[ind_gas]*m_gas[ind_gas]),
                             np.sum(l1_gas[ind_gas]*m_gas[ind_gas]),
                             np.sum(l2_gas[ind_gas]*m_gas[ind_gas]) ] )

         L_dir = L_gas / np.sqrt( np.dot(L_gas,L_gas) )
         j_gas = l0_gas[ind_gas] * L_dir[0] + l1_gas[ind_gas] * L_dir[1] + l2_gas[ind_gas] * L_dir[2]
         ind_bulge = np.where( j_gas < 0 )[0]
         MbulgeGas = np.sum( m_gas[ind_gas[ind_bulge]] ) * 2.

         for k in range(3):
            vel = np.sum( m_gas[ind_gas] * v_gas[ind_gas,k] ) / GasMass
            vel2 = np.sum( m_gas[ind_gas] * v_gas[ind_gas,k]**2 ) / GasMass
            SigmaGas[k] = np.sqrt( np.abs( vel2 - vel**2 ) )

         v_gas_z = v_gas[ind_gas,0] * L_dir[0] + v_gas[ind_gas,1] * L_dir[1] + v_gas[ind_gas,2] * L_dir[2]
         vel = np.sum( m_gas[ind_gas] * v_gas_z ) / GasMass
         vel2 = np.sum( m_gas[ind_gas] * v_gas_z**2 ) / GasMass
         SigmaGasZ = np.sqrt( np.abs( vel2 - vel**2 ) )

         phi_dir = np.cross( L_dir, r_gas, axisb=1 )
         phi_dir /= np.linalg.norm(phi_dir,axis=1)[:,np.newaxis]
         v_gas_phi = v_gas[ind_gas,0] * phi_dir[ind_gas,0] + v_gas[ind_gas,1] * phi_dir[ind_gas,1] + v_gas[ind_gas,2] * phi_dir[ind_gas,2]
         VgasPhi = np.sum( m_gas[ind_gas]*v_gas_phi ) / GasMass

         m_gas_perc[:] = np.percentile( m_gas[ind_gas], percentiles )

         ind_in = np.where( Vr_gas[ind_gas] < 0 )[0]
         if ind_in.size > 0:
           ind_in = ind_gas[ind_in]
           GasMassInflow = np.sum( m_gas[ind_in] )
           GasInflowVel = np.sum( Vr_gas[ind_in]*m_gas[ind_in] ) / GasMassInflow
           GasInflowVel_perc[:] = np.percentile( Vr_gas[ind_in], percentiles )

         ind_out = np.where( Vr_gas[ind_gas] >= 0 )[0]
         if ind_out.size > 0:
           ind_out = ind_gas[ind_out]
           GasMassOutflow = np.sum( m_gas[ind_out] )
           GasOutflowVel = np.sum( Vr_gas[ind_out]*m_gas[ind_out] ) / GasMassOutflow
           GasOutflowVel_perc[:] = np.percentile( Vr_gas[ind_out], percentiles )


       # --- save results ---

       var['prf']['Ndm'][n,i] = Ndm
       var['prf']['DarkMass'][n,i] = DarkMass

       var['prf']['Nstar'][n,i] = Nstar
       var['prf']['StarMass'][n,i] = StarMass
       var['prf']['MbulgeStar'][n,i] = MbulgeStar
       var['prf']['SigmaStar'][n,i,:] = SigmaStar
       var['prf']['L_star'][n,i,:] = L_star
       var['prf']['m_star_perc'][n,i,:] = m_star_perc

       var['prf']['Ngas'][n,i] = Ngas
       var['prf']['GasMass'][n,i] = GasMass
       var['prf']['MbulgeGas'][n,i] = MbulgeGas
       var['prf']['SFR'][n,i] = SFR
       var['prf']['SigmaGas'][n,i,:] = SigmaGas
       var['prf']['L_gas'][n,i,:] = L_gas
       var['prf']['m_gas_perc'][n,i,:] = m_gas_perc

       var['prf']['Temp'][n,i] = Temp
       var['prf']['SoundSpeed'][n,i] = SoundSpeed
       var['prf']['SigmaGasZ'][n,i] = SigmaGasZ
       var['prf']['VgasPhi'][n,i] = VgasPhi

       var['prf']['GasMassInflow'][n,i] = GasMassInflow
       var['prf']['GasInflowVel'][n,i] = GasInflowVel
       var['prf']['GasInflowVel_perc'][n,i,:] = GasInflowVel_perc
       var['prf']['GasMassOutflow'][n,i] = GasMassOutflow
       var['prf']['GasOutflowVel'][n,i] = GasOutflowVel
       var['prf']['GasOutflowVel_perc'][n,i,:] = GasOutflowVel_perc



    print snapnum, ':  redshift =', var['redshift'][n]
    sys.stdout.flush()


 ############################
  ### WRITE DATA TO FILE ###
 ############################

 #outfile = simdir + '/bhdata.hdf5'
 outfile = simdir + '/bh_rbins.hdf5'
 print '... writing', outfile
 if os.path.isfile(outfile ):
   os.remove(outfile)
 file = h5py.File(outfile, 'w')
 for grname in var.keys():
    if type(var[grname]) is dict:
      group = file.create_group( grname )
      for keyname in var[grname].keys():
         group.create_dataset(keyname, data=var[grname][keyname])
    else:
      file.create_dataset(grname, data=var[grname])
 file.close()





print '\n\nDone!!'

