# run as >>> execfile('bhfeedback_SFRhistory.py')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as cl
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
#import astropy.units as u
#from astropy.cosmology import FlatLambdaCDM, z_at_value

import daa_lib as daa
from daa_constants import *

import imp
imp.reload(daa)



#######################################
 ### DEFINE SIMULATION DIRECTORIES ###
#######################################

IC_list = [ 'sn155'#,
            #'sn214',
            #'sn277'
          ]

sim_list = [ 'nof_s8e1_n128'#,
             ##'spawn_s7e1v30f05_w100dt1000m100_TshVsh',
             #'spawn_s8e1v30f05_w100dt1000m100_TshVsh',
             #'spawn_s9e1v30f05_w100dt1000m100_TshVsh'
           ]


basedir = '/simons/scratch/dangles/FIRE/bhfeedback/h113_HR_'




Rgal = 10.  # kpc
dt_bin = 0.01  # Myr

for ICdir in IC_list:
   for sname in sim_list:

      simdir = basedir + ICdir + '/' + sname
      outdir = simdir + '/cenBH/'
      bh_file = outdir + 'centralBH.hdf5'   
      print '...reading ' + bh_file
      bh = daa.load_dict_from_hdf5( bh_file )
      snapnum = bh['snapnum'][-1]

      
      Pbh = g.readsnap( simdir, snapnum, 5, cosmological=1 )
      cen_pos = Pbh['p'][0,:]
      Pstar = g.readsnap( simdir, snapnum, 4, cosmological=1 )
      r_star = Pstar['p'][:,:] - cen_pos[:]
      R_star = np.sqrt((r_star*r_star).sum(axis=1))
      indR = np.where(R_star <= Rgal)[0]
      m_star = Pstar['m'][indR] * UnitMass_in_Msun
      aform = Pstar['age'][indR] 
      zform = 1./aform - 1.
      tform = 1e3 * util.age_of_universe( zform, h=Pbh['hubble'], Omega_M=Pbh['omega_matter'] )  # in Myr

      print 'snapnum =', snapnum, '  Mstar[Msun] =', m_star.sum()
      t_bins = np.arange( tform.min(), tform.max()+dt_bin, dt_bin )
      mass_in_bin, bin_edges, binum = binned_statistic( tform, m_star, statistic='sum', bins=t_bins )
      SFR = mass_in_bin / (dt_bin*1e6)
      t_avg = ( bin_edges[1:] + bin_edges[:-1] ) / 2.
      ind_recent = (t_avg >= bh['time'][0])

      ############################
       ### WRITE DATA TO FILE ###
      ############################
      var = { 'Rgal':Rgal, 'dt_bin':dt_bin, 'SFR':SFR, 'time':t_avg, 'SFR_recent':SFR[ind_recent], 'time_recent':t_avg[ind_recent]-t_avg[ind_recent].min() }
      outfile = outdir + 'SFH_%04d.hdf5' % snapnum
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
      

      # --- test SFR spike ...
      snapnum = bh['snapnum'][0]
      Pbh = g.readsnap( simdir, snapnum, 5, cosmological=1 )
      cen_pos = Pbh['p'][0,:]
      Pgas = g.readsnap( simdir, snapnum, 0, cosmological=1 )
      r_gas = Pgas['p'][:,:] - cen_pos[:]
      R_gas = np.sqrt((r_gas*r_gas).sum(axis=1))
      sfr = Pgas['sfr']
      rho = Pgas['rho'] * UnitDensity_in_cgs / PROTONMASS
      T = g.gas_temperature(Pgas['u'],Pgas['ne'])
      ind = np.where( (R_gas <= Rgal) & (sfr > 0.) )[0]      
      percentiles = [ 10, 50, 90, 99 ]

      sfr_tot = np.sum(sfr[ind])
      ind_sort = np.argsort( sfr[ind] )
      sfr_sorted = sfr[ind][ind_sort]
      rho_sorted = rho[ind][ind_sort]
      T_sorted = T[ind][ind_sort]
      R_sorted = R_gas[ind][ind_sort]

      sfr_cum = np.cumsum( sfr_sorted )
      ind_cum = np.where( sfr_cum > 0.1*sfr_tot )[0]
      
      print 'SFR = ', np.sum(sfr[ind]), '  perc =', np.percentile( sfr[ind], percentiles )


      Pgas_end = g.readsnap( simdir, bh['snapnum'][-1], 0, cosmological=1 )
      rho_end = Pgas_end['rho'] * UnitDensity_in_cgs / PROTONMASS
      sfr_end = Pgas_end['sfr']      
















print '\nDone!!\n'




