#>>> execfile('bhfeedback_gridTEST.py')

import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os as os
import sys as sys
import numpy as np
from scipy.stats import binned_statistic_2d
import h5py as h5py
import pandas as pd
from glob import glob
from copy import deepcopy
import gadget as g
import utilities as util
import visualization as vis
import daa_lib as daa
from daa_constants import *

import imp
imp.reload(daa)
imp.reload(util)



"""
narg = len(sys.argv)
if narg == 5:
  basedir = sys.argv[1]
  sname = sys.argv[2]
  snap_ini = int(sys.argv[3])
  snap_end = int(sys.argv[4])
else:
  print 'syntax: python bhfeedback_cenBHgrids.py basedir sname snap_ini snap_end'
  sys.exit()
"""


basedir = '/simons/scratch/dangles/FIRE/bhfeedback/'
sname = 'h113_HR_sn277/spawn_s9e1v30f05_w100dt1000m100_TshVsh'
snap_ini = 200
snap_end = 201


SkipIfDone = 1

npix = 512

dx0 = 7. 

bfr = 2.


simdir = basedir + sname

outdir = simdir + '/grids'
if not os.path.exists(outdir):
  try:
    os.mkdir(outdir)
  except OSError:
    try:
      os.makedirs(outdir)
    except OSError:
      print '...could not create output directory: ' + outdir


grid_vars = { 'mgas0':{ 'dx':dx0, 'Npart':0, 'g':np.zeros([npix,npix]) },
              'mgas_wind0':{ 'dx':dx0, 'Npart':0, 'g':np.zeros([npix,npix]) },
              'mgas_nowind0':{ 'dx':dx0, 'Npart':0, 'g':np.zeros([npix,npix]) },
              'T0':{ 'dx':dx0, 'Npart':0, 'g':np.zeros([npix,npix]) },
              'T_wind0':{ 'dx':dx0, 'Npart':0, 'g':np.zeros([npix,npix]) },
              'T_nowind0':{ 'dx':dx0, 'Npart':0, 'g':np.zeros([npix,npix]) }
            }

grid_tmp = { 'snapnum':-1, 'redshift':-1., 
             'faceon':deepcopy(grid_vars), 'edgeon':deepcopy(grid_vars) }



# --- rotation matrix ---
bhdata = daa.load_dict_from_hdf5(simdir + '/cenBH/cenBH_0001.hdf5')
#bhdata = daa.load_dict_from_hdf5(simdir + '/cenBH/cenBH_0032.hdf5')
Rref = 1.   # kpc
ind_Rcm = np.argmin( np.abs(bhdata['cum']['radius']-Rref) )
Lgas = bhdata['cum']['L_gas'][ind_Rcm,:]
zaxis = np.array([0.,0.,1.])
L_dir = Lgas / np.linalg.norm(Lgas)
rotax = np.cross(L_dir,zaxis)
rotangle = np.arccos(np.dot(L_dir,zaxis))
rotmatrix = daa.rotation_matrix( rotax, rotangle )


for snapnum in range(snap_ini,snap_end,1):
#for snapnum in [0,1]:

  outfile = outdir + '/gridTEST_%04d.hdf5' % snapnum
  if SkipIfDone and os.path.exists( outfile ):
    print "\nfile exists:  %s\n" % outfile
    continue

  try:
    header = g.readsnap( simdir, snapnum, 0, header_only=1)
  except Exception:
    print simdir, '  snapnum=', snapnum, ' not found!'
    continue

  if header['k'] == -1:
    print simdir, '  snapnum=', snapnum, ' not found!'
    continue

  grid = grid_tmp.copy()

  grid['snapnum'] = snapnum
  grid['redshift'] = header['redshift']


  # --- REFERENCE FRAME ---
  bhfile = simdir + '/cenBH/cenBH_%04d.hdf5' % snapnum
  if not os.path.exists( bhfile ):
    print ' --> snapshot ', snapnum, '  not available in bhdata! '
    continue
  bhdata = daa.load_dict_from_hdf5(bhfile)
  cen_pos = bhdata['cum']['COMstar_pos'][ind_Rcm,:] + bhdata['bh']['pos'][:]


  # --- GAS ---
  P = g.readsnap( simdir, snapnum, 0, cosmological=1 )
  ID_gas = P['id']
  r_gas = P['p'][:,:] - cen_pos[:]
  r_gas_rot = np.tensordot( rotmatrix, r_gas[:,:], axes=(1,1) ).T
  m = P['m'][:] * UnitMass_in_Msun
  T = g.gas_temperature(P['u'],P['ne'])
  hsml = P['h'][:]
  AGNWindID = 1913298393
  ind_wind = np.where( ID_gas == AGNWindID )[0]
  IsWind = ( m < 200. )
  IsNotWind = ( m > 200. )
  

  for im, dx in zip( ['0'], [dx0] ):
     xrg = [ -dx, dx ]

     # --- everything ---

     ind = np.where( (np.abs(r_gas_rot[:,0]) <= xrg[1]*bfr) & (np.abs(r_gas_rot[:,1]) <= xrg[1]*bfr) & (np.abs(r_gas_rot[:,2]) <= xrg[1]*bfr) )[0]
     if ind.size <= 1:
       continue
     x_grid = r_gas_rot[ind,0]
     y_grid = r_gas_rot[ind,1]
     z_grid = r_gas_rot[ind,2]
     m_grid = m[ind]
     T_grid = T[ind]
     hsml_grid = daa.get_particle_hsml( x_grid, y_grid, z_grid, DesNgb=32 )

     grid['faceon']['mgas'+im]['Npart'] = ind.size
     grid['faceon']['T'+im]['Npart'] = ind.size
     m_sMap, Tm_sMap = daa.make_2Dgrid( x_grid, y_grid, xrg, weight1=m_grid, weight2=T_grid*m_grid, hsml=hsml_grid, pixels=npix )
     Mgas_map, Tgas_map = daa.clip_2Dgrid( m_sMap, map2=Tm_sMap )
     grid['faceon']['mgas'+im]['g'][:,:] = Mgas_map.T
     grid['faceon']['T'+im]['g'][:,:] = Tgas_map.T

     grid['edgeon']['mgas'+im]['Npart'] = ind.size
     grid['edgeon']['T'+im]['Npart'] = ind.size
     m_sMap, Tm_sMap = daa.make_2Dgrid( x_grid, z_grid, xrg, weight1=m_grid, weight2=T_grid*m_grid, hsml=hsml_grid, pixels=npix )
     Mgas_map, Tgas_map = daa.clip_2Dgrid( m_sMap, map2=Tm_sMap )
     grid['edgeon']['mgas'+im]['g'][:,:] = Mgas_map.T
     grid['edgeon']['T'+im]['g'][:,:] = Tgas_map.T

     # --- normal gas (not wind) ---

     ind = np.where( IsNotWind & (np.abs(r_gas_rot[:,0]) <= xrg[1]*bfr) & (np.abs(r_gas_rot[:,1]) <= xrg[1]*bfr) & (np.abs(r_gas_rot[:,2]) <= xrg[1]*bfr) )[0]
     if ind.size <= 1:
       continue
     x_grid = r_gas_rot[ind,0]
     y_grid = r_gas_rot[ind,1]
     z_grid = r_gas_rot[ind,2]
     m_grid = m[ind]
     T_grid = T[ind]
     hsml_grid = daa.get_particle_hsml( x_grid, y_grid, z_grid, DesNgb=32 )

     grid['faceon']['mgas_nowind'+im]['Npart'] = ind.size
     grid['faceon']['T_nowind'+im]['Npart'] = ind.size
     m_sMap, Tm_sMap = daa.make_2Dgrid( x_grid, y_grid, xrg, weight1=m_grid, weight2=T_grid*m_grid, hsml=hsml_grid, pixels=npix )
     Mgas_map, Tgas_map = daa.clip_2Dgrid( m_sMap, map2=Tm_sMap )
     grid['faceon']['mgas_nowind'+im]['g'][:,:] = Mgas_map.T
     grid['faceon']['T_nowind'+im]['g'][:,:] = Tgas_map.T

     grid['edgeon']['mgas_nowind'+im]['Npart'] = ind.size
     grid['edgeon']['T_nowind'+im]['Npart'] = ind.size
     m_sMap, Tm_sMap = daa.make_2Dgrid( x_grid, z_grid, xrg, weight1=m_grid, weight2=T_grid*m_grid, hsml=hsml_grid, pixels=npix )
     Mgas_map, Tgas_map = daa.clip_2Dgrid( m_sMap, map2=Tm_sMap )
     grid['edgeon']['mgas_nowind'+im]['g'][:,:] = Mgas_map.T
     grid['edgeon']['T_nowind'+im]['g'][:,:] = Tgas_map.T

     # --- wind particles ---

     ind = np.where( IsWind & (np.abs(r_gas_rot[:,0]) <= xrg[1]*bfr) & (np.abs(r_gas_rot[:,1]) <= xrg[1]*bfr) & (np.abs(r_gas_rot[:,2]) <= xrg[1]*bfr) )[0]
     if ind.size <= 1:
       continue
     x_grid = r_gas_rot[ind,0]
     y_grid = r_gas_rot[ind,1]
     z_grid = r_gas_rot[ind,2]
     m_grid = m[ind]
     T_grid = T[ind]
     hsml_grid = daa.get_particle_hsml( x_grid, y_grid, z_grid, DesNgb=32 )

     grid['faceon']['mgas_wind'+im]['Npart'] = ind.size
     grid['faceon']['T_wind'+im]['Npart'] = ind.size
     m_sMap, Tm_sMap = daa.make_2Dgrid( x_grid, y_grid, xrg, weight1=m_grid, weight2=T_grid*m_grid, hsml=hsml_grid, pixels=npix )
     Mgas_map, Tgas_map = daa.clip_2Dgrid( m_sMap, map2=Tm_sMap )
     grid['faceon']['mgas_wind'+im]['g'][:,:] = Mgas_map.T
     grid['faceon']['T_wind'+im]['g'][:,:] = Tgas_map.T

     grid['edgeon']['mgas_wind'+im]['Npart'] = ind.size
     grid['edgeon']['T_wind'+im]['Npart'] = ind.size
     m_sMap, Tm_sMap = daa.make_2Dgrid( x_grid, z_grid, xrg, weight1=m_grid, weight2=T_grid*m_grid, hsml=hsml_grid, pixels=npix )
     Mgas_map, Tgas_map = daa.clip_2Dgrid( m_sMap, map2=Tm_sMap )
     grid['edgeon']['mgas_wind'+im]['g'][:,:] = Mgas_map.T
     grid['edgeon']['T_wind'+im]['g'][:,:] = Tgas_map.T


  # --- STARS ---  
  """
  P = g.readsnap( simdir, snapnum, 4, cosmological=1 )
  x = P['p'][:,0] - cen_pos[0]
  y = P['p'][:,1] - cen_pos[1]
  z = P['p'][:,2] - cen_pos[2]
  m = P['m'][:] * UnitMass_in_Msun
  hsml = daa.get_particle_hsml( x, y, z, DesNgb=32 )
  for im, dx in zip( ['1','2','3'], [dx1,dx2,dx3] ):
     xrg = [ -dx, dx ]
     ind = np.where( (np.abs(x) <= xrg[1]*bfr) & (np.abs(y) <= xrg[1]*bfr) & (np.abs(z) <= xrg[1]*bfr) )[0]
     grid['mstar'+im]['Npart'] = ind.size
     if ind.size <= 1:
       continue
     x_grid = x[ind]
     y_grid = y[ind]
     m_grid = m[ind]
     hsml_grid = hsml[ind]
     m_sMap = daa.make_2Dgrid( x_grid, y_grid, xrg, weight1=m_grid, hsml=hsml_grid, pixels=npix )
     Mstar_map = daa.clip_2Dgrid( m_sMap )
     grid['mstar'+im]['g'][:,:] = Mstar_map.T
  """

  ############################
   ### WRITE DATA TO FILE ###
  ############################

  print '\n... writing', outfile
  if os.path.isfile(outfile ):
    os.remove(outfile)
  file = h5py.File(outfile, 'w')
  for grname in grid.keys():
    if type(grid[grname]) is dict:
      group = file.create_group( grname )
      for keyname in grid[grname].keys():
         if type(grid[grname][keyname]) is dict:
           group2 = group.create_group( keyname )
           for keyname2 in grid[grname][keyname].keys():
              group2.create_dataset(keyname2, data=grid[grname][keyname][keyname2])
         else:
           group.create_dataset(keyname, data=grid[grname][keyname])
    else:
      file.create_dataset(grname, data=grid[grname])
  file.close()



###################
 ### MAKE PLOT ###
###################

outdir = './bhfeedback/' + sname + '/gridTEST'
if not os.path.exists(outdir):
  try:
    os.mkdir(outdir)
  except OSError:
    try:
      os.makedirs(outdir)
    except OSError:
      print '...could not create output directory: ' + outdir

for snapnum in range(snap_ini,snap_end,1):

  outfile = outdir + '/gridTEST_%04d.png' % snapnum
  gridfile = simdir + '/grids/gridTEST_%04d.hdf5' % snapnum
  if not os.path.exists( gridfile ):
    print "\nfile not found:  %s\n" % gridfile
    continue
  grd = daa.load_dict_from_hdf5( gridfile )

  grd = grd['edgeon']

  print '\n', gridfile, ' --> ', outfile
  fig = plt.figure( figsize=(20.,5.), dpi=150 )
  fig.subplots_adjust(left=0., right=1., bottom=0., top=1., hspace=0., wspace=0.)
  ax1 = fig.add_subplot(141)
  ax2 = fig.add_subplot(142)
  ax3 = fig.add_subplot(143)
  ax4 = fig.add_subplot(144)

  vmin = [ 6.5, 6.5, 3., 6.5 ]
  vmax = [ 9., 9., 7., 9. ]
  imgs = ['mgas0', 'mgas_nowind0', 'mgas_wind0', 'mgas_nowind0']
  for j, ax in enumerate([ax1,ax2,ax3,ax4]):
    img_str = imgs[j]
    dx = grd[img_str]['dx']
    xrg = [ -dx, dx ]
    extent = [xrg[0],xrg[1],xrg[0],xrg[1]]
    ax.set_xlim(xrg); ax.set_ylim(xrg);
    ax.set_aspect('equal'); ax.axis('off')
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
       item.set_visible(False)
    if j == 3:
      im_temp = ax.imshow( np.log10( grd['mgas_nowind0']['g'][:,:]+grd['mgas_wind0']['g'][:,:] ), vmin=vmin[j], vmax=vmax[j], cmap=cm.Spectral, interpolation='bicubic', origin='low', extent=extent )
    else:
      im_temp = ax.imshow( np.log10( grd[img_str]['g'][:,:] ), vmin=vmin[j], vmax=vmax[j], cmap=cm.Spectral, interpolation='bicubic', origin='low', extent=extent )

  fig.savefig( outfile, dpi=150 )
  plt.close()






print '\nDone!\n'




