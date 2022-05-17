# >>> execfile('bhfeedback_cenBH_combine.py')

import numpy as np
import gadget as g
import os as os
import sys as sys
import h5py as h5py
from glob import glob
from copy import deepcopy
import daa_lib as daa



#####################################################
def data_template(nsnap,nbin,nperc):
 var = {}
 var['radial_bins'] = np.zeros(nbin+1)
 var['snapnum'] = np.zeros(nsnap,dtype=np.int32)
 var['redshift'] = np.zeros(nsnap)
 var['time'] = np.zeros(nsnap)
 var['cum'] = { 'Vcirc':np.zeros([nsnap,nbin]),
                'Ndm':np.zeros([nsnap,nbin],dtype=np.int32), 'DarkMass':np.zeros([nsnap,nbin]),
                'Nstar':np.zeros([nsnap,nbin],dtype=np.int32), 'StarMass':np.zeros([nsnap,nbin]), 'L_star':np.zeros([nsnap,nbin,3]), 'MbulgeStar':np.zeros([nsnap,nbin]),
                'COMstar_pos':np.zeros([nsnap,nbin,3]), 'COMstar_vel':np.zeros([nsnap,nbin,3]), 'SigmaStar':np.zeros([nsnap,nbin,3]),
                'Ngas':np.zeros([nsnap,nbin],dtype=np.int32), 'GasMass':np.zeros([nsnap,nbin]), 'L_gas':np.zeros([nsnap,nbin,3]), 'MbulgeGas':np.zeros([nsnap,nbin]),
                'COMgas_pos':np.zeros([nsnap,nbin,3]), 'COMgas_vel':np.zeros([nsnap,nbin,3]), 'SigmaGas':np.zeros([nsnap,nbin,3]),
                'SFR':np.zeros([nsnap,nbin]),
                'Ngas_BAL':np.zeros([nsnap,nbin],dtype=np.int32), 'GasMass_BAL':np.zeros([nsnap,nbin]), 'L_gas_BAL':np.zeros([nsnap,nbin,3]), 'COMgas_pos_BAL':np.zeros([nsnap,nbin,3]), 'COMgas_vel_BAL':np.zeros([nsnap,nbin,3])
              }
 var['sph'] = { 'Ndm':np.zeros([nsnap,nbin],dtype=np.int32), 'DarkMass':np.zeros([nsnap,nbin]), 'm_dm_perc':np.zeros([nsnap,nbin,nperc]),
                'Nstar':np.zeros([nsnap,nbin],dtype=np.int32), 'StarMass':np.zeros([nsnap,nbin]), 'SigmaStar':np.zeros([nsnap,nbin,3]), 'VstarPhi':np.zeros([nsnap,nbin]),
                'L_star':np.zeros([nsnap,nbin,3]), 'm_star_perc':np.zeros([nsnap,nbin,nperc]),
                'Ngas':np.zeros([nsnap,nbin],dtype=np.int32), 'GasMass':np.zeros([nsnap,nbin]), 'SigmaGas':np.zeros([nsnap,nbin,3]),
                'L_gas':np.zeros([nsnap,nbin,3]), 'm_gas_perc':np.zeros([nsnap,nbin,nperc]),
                'Temp':np.zeros([nsnap,nbin]), 'SoundSpeed':np.zeros([nsnap,nbin]), 'VgasPhi':np.zeros([nsnap,nbin]),
                'GasMassInflow':np.zeros([nsnap,nbin]), 'GasInflowVel':np.zeros([nsnap,nbin]), 'GasInflowVel_perc':np.zeros([nsnap,nbin,nperc]), 'GasInflowRate':np.zeros([nsnap,nbin]),
                'GasMassOutflow':np.zeros([nsnap,nbin]), 'GasOutflowVel':np.zeros([nsnap,nbin]), 'GasOutflowVel_perc':np.zeros([nsnap,nbin,nperc]), 'GasOutflowRate':np.zeros([nsnap,nbin]),
                'SFR':np.zeros([nsnap,nbin]),
                'Ngas_BAL':np.zeros([nsnap,nbin],dtype=np.int32), 'GasMass_BAL':np.zeros([nsnap,nbin]), 'm_gas_perc_BAL':np.zeros([nsnap,nbin,nperc]), 'L_gas_BAL':np.zeros([nsnap,nbin,3]), 'Temp_BAL':np.zeros([nsnap,nbin]), 'SoundSpeed_BAL':np.zeros([nsnap,nbin]),
                'GasMassOutflow_BAL':np.zeros([nsnap,nbin]), 'GasOutflowVel_BAL':np.zeros([nsnap,nbin]), 'GasOutflowVel_perc_BAL':np.zeros([nsnap,nbin,nperc]), 'GasOutflowRate_BAL':np.zeros([nsnap,nbin])
              }
 cylstr = { 'Nstar':np.zeros([nsnap,nbin],dtype=np.int32), 'StarMass':np.zeros([nsnap,nbin]), 'SigmaStar':np.zeros([nsnap,nbin,3]), 'L_star':np.zeros([nsnap,nbin,3]), 'VstarPhi':np.zeros([nsnap,nbin]),
            'Ngas':np.zeros([nsnap,nbin],dtype=np.int32), 'GasMass':np.zeros([nsnap,nbin]), 'SigmaGas':np.zeros([nsnap,nbin,3]), 'L_gas':np.zeros([nsnap,nbin,3]),
            'Temp':np.zeros([nsnap,nbin]), 'SoundSpeed':np.zeros([nsnap,nbin]), 'VgasPhi':np.zeros([nsnap,nbin]),
            'GasMassInflow':np.zeros([nsnap,nbin]), 'GasInflowVel':np.zeros([nsnap,nbin]), 'GasInflowVel_perc':np.zeros([nsnap,nbin,nperc]), 'GasInflowRate':np.zeros([nsnap,nbin]),
            'GasMassOutflow':np.zeros([nsnap,nbin]), 'GasOutflowVel':np.zeros([nsnap,nbin]), 'GasOutflowVel_perc':np.zeros([nsnap,nbin,nperc]), 'GasOutflowRate':np.zeros([nsnap,nbin]),
            'SFR':np.zeros([nsnap,nbin])
           }
 z_cyl = [ '100pc', '1kpc' ]
 var['cyl'] = {}
 for ztag in z_cyl:
    var['cyl'][ztag] = deepcopy(cylstr)
 var['bh'] = { 'Mass':np.zeros(nsnap), 'BH_Mass':np.zeros(nsnap), 'pos':np.zeros([nsnap,3]), 'vel':np.zeros([nsnap,3]), 'RadiusOfInfluence':np.zeros(nsnap) }
 return var
#####################################################



IC_list = [ #'h113_HR_sn155',
            #'h113_HR_sn214',
            #'h113_HR_sn277',
            'h113_HR_sn152'
          ]

sim_list = [ #'spawn_s9e1v30f05_w100dt1000m100_TshVsh',
             #'spawn_s8e1v30f05_w100dt1000m100_TshVsh',
             ##'spawn_s7e1v30f05_w100dt1000m100_TshVsh',
             #'nof_s8e1_n128'#,
             #'spawn_s9e1v30f05_w100dt1000m100_TshVsh_MERGE',
             #'spawn_s9e1v30f05_w10dt100m1000_TshVsh',
             #'spawn_s9e1v30f05_w100dt1000m100',
             #'spawn_s9e1v30f033_w100dt1000m100_TshVsh',
             #'spawn_s9e1v30f009_w100dt1000m100_TshVsh',
             #'kick_s9e1v15000p2_n128#',
             #'kick_s9e1v1500p20_n128',
             #'bal_s9e1v30f05_n128',
             #'TEST_SnapOrig_nof_s8e1_n128',
             #'TEST_GizmoOrig_nof_s8e1_n128',
             #'TEST_noMetDiff_nof_s8e1_n128'
             #'kick_COL_s9e1v15000p2_n128',
             #'kick_COL_s9e1v1500p20_n128'
#             'nof_s8e1_n128',
#             'spawn_res251_s9e01v30f05_w10dt100m1000_TshVsh',
#             'spawn_res251_s9e1v30f05_w10dt100m1000_TshVsh',
#             'spawn_res251_s9e1v30f033_w10dt100m1000_TshVsh',
#             'spawn_res251_s9e1v30f02_w10dt100m1000_TshVsh',
#             'spawn_res251_s9e1v30f011_w10dt100m1000_TshVsh',
#             'spawn_res251_s9e1v30f009_w10dt100m1000_TshVsh'
#             'TEST_mainbranch_e1f05',
#             'TEST_mainbranch_e1f009'
             'spawn_res251_s9e01v30f05_w10dt100m1000_TshVsh',
             'spawn_res251_s9e1v30f05_w10dt100m1000_TshVsh',
             'spawn_res251_s9e1v30f033_w10dt100m1000_TshVsh',
             'spawn_res251_s9e1v30f02_w10dt100m1000_TshVsh',
             'spawn_res251_s9e1v30f009_w10dt100m1000_TshVsh'
           ]

basedir = '/simons/scratch/dangles/FIRE/bhfeedback/'


#isim = 15
#sim_list = [ sim_list[isim] ]

nsim = len(sim_list)

for ICdir in IC_list:
 for sname in sim_list:

  simdir = basedir + ICdir + '/' + sname
  outdir = simdir + '/cenBH'

  # --- read first cenBH file and define output structure:
  snaplist, filestr = daa.file_list( outdir + '/cenBH', ext='.hdf5' )
  nsnap = snaplist.size
  print snaplist 

  cenBH_1 = daa.load_dict_from_hdf5(outdir + '/cenBH_%04d.hdf5' % snaplist[0])
  nbin, nperc = cenBH_1['sph']['GasInflowVel_perc'].shape
  var = data_template(nsnap,nbin,nperc)
  var['radial_bins'][:] = cenBH_1['sph']['radial_bins'][:]
  
  for n in range(nsnap):

     snapnum = snaplist[n]
     cenBH_path = outdir + '/cenBH_%04d.hdf5' % snapnum
     if not os.path.exists( cenBH_path ):
       print "\nfile doesn't exist!!  %s\n" % cenBH_path
       break

     cenBH = daa.load_dict_from_hdf5(cenBH_path)

     for key1 in var.keys():
       if key1 in ['radial_bins']:
         continue 

       if type(var[key1]) is dict:
      
         for key2 in var[key1].keys():

            if type(var[key1][key2]) is dict:

              for key3 in var[key1][key2].keys():

                 if type(var[key1][key2][key3]) is dict:
                   print '....what is this?', key1, key2, key3
                 else:
                   if len(var[key1][key2][key3].shape) == 1:
                     var[key1][key2][key3][n] = cenBH[key1][key2][key3]
                   else:
                     var[key1][key2][key3][n,:] = cenBH[key1][key2][key3][:]

            else:
              if len(var[key1][key2].shape) == 1:
                var[key1][key2][n] = cenBH[key1][key2]
              else:
                var[key1][key2][n,:] = cenBH[key1][key2][:]

       else:
         var[key1][n] = cenBH[key1]


  ############################
   ### WRITE DATA TO FILE ###
  ############################

  outfile = outdir + '/centralBH.hdf5'
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




print '\n\nDone!\n'






