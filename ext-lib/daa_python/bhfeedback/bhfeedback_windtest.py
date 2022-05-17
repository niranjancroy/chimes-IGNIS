# run as >>> execfile('bhfeedback_windtest.py')
# python bhfeedback_windtest.py  > log/windtest.log 2>&1 &

import numpy as np
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

np.set_printoptions(linewidth=150)


simdir = '/simons/scratch/dangles/FIRE/bhfeedback/h113_HR_sn152/spawn_res251_s9e1v30f05_w10dt100m1000_TshVsh'
snaplist = range(260,370)
ID_BAL = 4060782041

for snapnum in snaplist:

 for simdir in [#'/simons/scratch/dangles/FIRE/bhfeedback/h113_HR_sn277/spawn_s9e1v30f05_w100dt1000m100_TshVsh',
                #'/simons/scratch/dangles/FIRE/bhfeedback/h113_HR_sn277/spawn_s9e1v30f05_w10dt100m1000_TshVsh',
                #'/simons/scratch/dangles/FIRE/bhfeedback/h113_HR_sn277/spawn_s9e1v30f05_w100dt1000m100'#,
                #'/simons/scratch/dangles/FIRE/bhfeedback/h113_HR_sn277/spawn_s9e1v30f05_w100dt1000m100_TshVsh_MERGE'#,
                #'/simons/scratch/dangles/FIRE/bhfeedback/h113_HR_sn152/nof_s8e1_n128'
                #'/simons/scratch/dangles/FIRE/bhfeedback/h113_HR_sn152/spawn_res251_s9e1v30f05_w10dt100m1000_TshVsh',
                #'/simons/scratch/dangles/FIRE/bhfeedback/h113_HR_sn152/TEST_mainbranch_e1f05'
                #'/simons/scratch/dangles/FIRE/bhfeedback/h113_HR_sn152/TEST_mainbranch_e1f009',
                '/simons/scratch/dangles/FIRE/bhfeedback/h113_HR_sn152/spawn_res251_s9e1v30f009_w10dt100m1000_TshVsh',
                '/simons/scratch/dangles/FIRE/bhfeedback/h113_HR_sn152/TEST_KungYi_e1f009' 
                #'/simons/scratch/dangles/FIRE/bhfeedback/h113_HR_sn277/kick_COL_s9e1v15000p2_n128'
               ]:

  print '\n', simdir, '     snapnum =', snapnum

  Pbh = g.readsnap( simdir, snapnum, 5, cosmological=1 )
  cen_pos = Pbh['p'][0,:]
  cen_vel = Pbh['v'][0,:]
  h = Pbh['hubble']
  hubble_factor = daa.hubble_z( Pbh['redshift'] )

  Pgas = g.readsnap( simdir, snapnum, 0, cosmological=1 )
  id_gas = Pgas['id']

  print '--- regular GAS ---'
  ind = np.where( id_gas != ID_BAL )[0]

  r_gas = Pgas['p'][ind,:] - cen_pos[:]
  R_gas = np.sqrt((r_gas*r_gas).sum(axis=1))

  m_gas = Pgas['m'][ind] * UnitMass_in_Msun
  v_gas = Pgas['v'][ind,:] - cen_vel[:] +  hubble_factor * r_gas * UnitLength_in_cm  / UnitVelocity_in_cm_per_s
  Vr_gas = (v_gas*r_gas).sum(axis=1) / R_gas
  T_gas = g.gas_temperature(Pgas['u'][ind],Pgas['ne'][ind])
  cs_gas = np.sqrt( GAMMA*GAMMA_MINUS1 * Pgas['u'][ind] )

  print 'm_gas =', m_gas.min(), m_gas.max()
  print 'Vr_gas =', Vr_gas.min(), Vr_gas.max()
  print 'T_gas =', T_gas.min(), T_gas.max()
  print 'R_gas =', R_gas.min(), R_gas.max()


  print '--- BAL winds ---'
  ind = np.where( id_gas == ID_BAL )[0]
  if ind.size == 0:
    continue

  r_wind = Pgas['p'][ind,:] - cen_pos[:]
  R_wind = np.sqrt((r_wind*r_wind).sum(axis=1))

  m_wind = Pgas['m'][ind] * UnitMass_in_Msun
  v_wind = Pgas['v'][ind,:] - cen_vel[:] +  hubble_factor * r_wind * UnitLength_in_cm  / UnitVelocity_in_cm_per_s
  Vr_wind = (v_wind*r_wind).sum(axis=1) / R_wind
  T_wind = g.gas_temperature(Pgas['u'][ind],Pgas['ne'][ind])
  cs_wind = np.sqrt( GAMMA*GAMMA_MINUS1 * Pgas['u'][ind] )

  print 'm_wind =', m_wind.min(), m_wind.max()
  print 'Vr_wind =', Vr_wind.min(), Vr_wind.max()
  print 'T_wind =', T_wind.min(), T_wind.max()
  print 'R_wind =', R_wind.min(), R_wind.max()

  if Vr_wind.max() > 5e4:
    print '********************************************************* look at this wind! *******************************************************************'

  ind_vel = np.where(np.abs(Vr_wind) > 0.99*Vr_wind.max())[0]
  if ind_vel.size == 0:
    continue
  print 'Vr_wind[ind_vel] =', Vr_wind[ind_vel]
  print 'R_wind[ind_vel] =', R_wind[ind_vel]
  print 'm_wind[ind_vel] =', m_wind[ind_vel]
  print 'T_wind[ind_vel] =', T_wind[ind_vel]
  #print 'cs_wind[ind_vel] =', cs_wind[ind_vel]

  sys.stdout.flush()





print '\n\nDone!\n'









