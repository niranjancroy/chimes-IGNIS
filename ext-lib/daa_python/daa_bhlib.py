
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#from matplotlib import rc
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



def bhrotmatrix( bhdata, Rref=1., inum=0 ):
   ind_Rcm = np.argmin( np.abs(bhdata['radial_bins']-Rref) )
   Lgas = bhdata['cum']['L_gas'][inum,ind_Rcm,:]
   zaxis = np.array([0.,0.,1.])
   L_dir = Lgas / np.linalg.norm(Lgas)
   rotax = np.cross(L_dir,zaxis)
   rotangle = np.arccos(np.dot(L_dir,zaxis))
   rotmatrix = daa.rotation_matrix( rotax, rotangle )
   return rotmatrix



def make_windgrid( simdir, snapnum, outfile, dx=1., rotmatrix=np.identity(3), center='BH', m_wind=3000., DesNgb=32, npix=512, SkipIfDone=0, P=-1 ):

   if SkipIfDone and os.path.exists( outfile ):
    print "\nfile exists:  %s\n" % outfile
    return 1

   outdir = simdir + '/tzgrid'
   if not os.path.exists(outdir):
     try:
       os.mkdir(outdir)
     except OSError:
       try:
         os.makedirs(outdir)
       except OSError:
         print '...could not create output directory: ' + outdir

   # --- REFERENCE FRAME ---
   bhfile = simdir + '/cenBH/cenBH_%04d.hdf5' % snapnum
   if not os.path.exists( bhfile ):
     print ' --> snapshot ', snapnum, '  not available in bhdata! '
     return 0
   bhdata = daa.load_dict_from_hdf5(bhfile)
   # --- BH as reference frame
   cen_pos = bhdata['bh']['pos'][:]
   cen_vel = bhdata['bh']['vel'][:]
   if center == 'COM':
     # --- COM of the stars
     Rref = 1. # kpc
     ind_Rcm = np.argmin( np.abs(bhdata['cum']['radius']-Rref) )
     cen_pos += bhdata['cum']['COMstar_pos'][ind_Rcm,:]  
     cen_vel += bhdata['cum']['COMstar_vel'][ind_Rcm,:] 

   grid_vars = { 'snapnum':snapnum, 'redshift':bhdata['redshift'], 'time':bhdata['time'], 'dx':dx, 'rotmatrix':rotmatrix,
                 'gas':{ 'Npart':0, 'm':np.zeros([npix,npix]), 'T':np.zeros([npix,npix]), 'sfr':np.zeros([npix,npix]) },
                 'wind':{ 'Npart':0, 'm':np.zeros([npix,npix]), 'T':np.zeros([npix,npix]), 'Vr':np.zeros([npix,npix]), 'vz':np.zeros([npix,npix]), 'ek':np.zeros([npix,npix]) }
               }
   grid = { 'faceon':deepcopy(grid_vars), 'edgeon':deepcopy(grid_vars) }

   if type(P) == 'dict':
     print ' ...making grid with given P:', simdir, 'snapnum, dx =', snapnum, dx
   else:
     print ' ...making grid:', simdir, '  [ snapnum, dx =', snapnum, dx, ']'
     P = g.readsnap( simdir, snapnum, 0, cosmological=1 )
   r_gas = P['p'][:,:] - cen_pos[:]
   R_gas = np.sqrt((r_gas*r_gas).sum(axis=1))
   hubble_factor = daa.hubble_z( P['redshift'] )
   v_gas = P['v'][:,:] - cen_vel[:] +  hubble_factor * r_gas * UnitLength_in_cm  / UnitVelocity_in_cm_per_s
   Vr_gas = (v_gas*r_gas).sum(axis=1) / R_gas
   m = P['m'][:] * UnitMass_in_Msun
   T = g.gas_temperature(P['u'],P['ne'])
   sfr = P['sfr'][:]
   hsml = P['h'][:]
   r_gas_rot = np.tensordot( rotmatrix, r_gas[:,:], axes=(1,1) ).T
   v_gas_rot = np.tensordot( rotmatrix, v_gas[:,:], axes=(1,1) ).T

   xrg = [ -dx, dx ]
   bfr = 2.

   # --- normal gas (not wind) ---
   ind = np.where( ( m > m_wind ) & (np.abs(r_gas_rot[:,0]) <= xrg[1]*bfr) & (np.abs(r_gas_rot[:,1]) <= xrg[1]*bfr) & (np.abs(r_gas_rot[:,2]) <= xrg[1]*bfr) )[0]
   if ind.size >= 1:
    x_grid = r_gas_rot[ind,0]
    y_grid = r_gas_rot[ind,1]
    z_grid = r_gas_rot[ind,2]
    m_grid = m[ind]
    T_grid = T[ind]
    sfr_grid = sfr[ind]
    hsml_grid = daa.get_particle_hsml( x_grid, y_grid, z_grid, DesNgb=DesNgb )
    for los in ['faceon','edgeon']:
       grid[los]['gas']['Npart'] = ind.size
       if los == 'faceon':
         coord_grid = y_grid
       else:
         coord_grid = z_grid
       m_sMap, Tm_sMap, sfr_sMap = daa.make_2Dgrid( x_grid, coord_grid, xrg, weight1=m_grid, weight2=T_grid*m_grid, weight3=sfr_grid, hsml=hsml_grid, pixels=npix )
       Mgas_map, Tgas_map = daa.clip_2Dgrid( m_sMap, map2=Tm_sMap )
       grid[los]['gas']['m'][:,:] = Mgas_map.T
       grid[los]['gas']['T'][:,:] = Tgas_map.T
       grid[los]['gas']['sfr'][:,:] = sfr_sMap.T

   # --- wind particles ---
   ind = np.where( ( m <= m_wind ) & (np.abs(r_gas_rot[:,0]) <= xrg[1]*bfr) & (np.abs(r_gas_rot[:,1]) <= xrg[1]*bfr) & (np.abs(r_gas_rot[:,2]) <= xrg[1]*bfr) )[0]
   if ind.size >= 1:
    x_grid = r_gas_rot[ind,0]
    y_grid = r_gas_rot[ind,1]
    z_grid = r_gas_rot[ind,2]
    Vr_grid = Vr_gas[ind]
    m_grid = m[ind]
    T_grid = T[ind]
    hsml_grid = daa.get_particle_hsml( x_grid, y_grid, z_grid, DesNgb=DesNgb/2 )
    for los in ['faceon','edgeon']:
       grid[los]['wind']['Npart'] = ind.size
       if los == 'faceon':
         coord_grid = y_grid
       else:
         coord_grid = z_grid
       m_sMap, Tm_sMap, Vr_sMap = daa.make_2Dgrid( x_grid, coord_grid, xrg, weight1=m_grid, weight2=T_grid*m_grid, weight3=Vr_grid*m_grid, hsml=hsml_grid, pixels=npix )
       Mgas_map, Tgas_map = daa.clip_2Dgrid( m_sMap, map2=Tm_sMap )
       Mgas_map, Vr_map  = daa.clip_2Dgrid( m_sMap, map2=Vr_sMap )
       grid[los]['wind']['m'][:,:] = Mgas_map.T
       grid[los]['wind']['T'][:,:] = Tgas_map.T
       grid[los]['wind']['Vr'][:,:] = Vr_map.T


   ############################
    ### WRITE DATA TO FILE ###
   ############################

   print '\n... writing', outfile
   if os.path.isfile(outfile):
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
   return 1


def make_windgrid_plot( gridfile, outfile, los='edgeon', vmax=-999 ):
   try:
     thegrid = daa.load_dict_from_hdf5( gridfile )
     print '\n', gridfile, ' --> ', outfile, '\n'
   except IOError:
     print 'grid not found: ', gridfile
     return 0, vmax
   grd = thegrid[los]

   fig = plt.figure( figsize=(5.,5.), dpi=150 )
   fig.subplots_adjust(left=0., right=1., bottom=0., top=1., hspace=0., wspace=0.)
   ax = fig.add_subplot(111)
   dx = grd['dx']
   xrg = [ -dx, dx ]
   extent = [xrg[0],xrg[1],xrg[0],xrg[1]]
   ax.set_xlim(xrg); ax.set_ylim(xrg);
   ax.set_aspect('equal'); ax.axis('off')
   for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
      item.set_visible(False)
   cmap_gas = cm.Greys
   cmap_wind = cm.jet
   alpha_wind = 0.25

   # --- background gas density
   gas_map = grd['gas']['m'][:,:]
   vmaximum = np.log10( gas_map.max() ) - 0.5
   if vmax == -999:
     vmax = vmaximum
   vmin = vmax - 3
   vgas = [ vmin, vmax ]
   #vgas = [ 5.5, 9.5 ]
   im = ax.imshow( np.log10( gas_map ), vmin=vgas[0], vmax=vgas[1], cmap=cmap_gas, interpolation='bicubic', origin='low', extent=extent )
   # --- AGN wind
   vwind = [ 3.5, 6 ]
   wind_grid = grd['wind']['m'][:,:]
   wind_map = np.ma.masked_array( wind_grid, wind_grid <= wind_grid.min() )
   im = ax.imshow( np.log10( wind_map ), vmin=vwind[0], vmax=vwind[1], cmap=cmap_wind, interpolation='bicubic', origin='low', extent=extent, alpha=alpha_wind )
   # --- dense ISM
   dense_map = np.ma.masked_array( gas_map, np.log10(gas_map) <= 8.5 )
   im = ax.imshow( np.log10( dense_map ), vmin=vgas[0], vmax=vgas[1], cmap=cmap_gas, interpolation='bicubic', origin='low', extent=extent, alpha=0.15 )
   # --- virial radius...
   acirc = np.linspace(0,2*np.pi,100)
   xcirc = np.cos(acirc)
   ycirc = np.sin(acirc)
   for rfac in [1.,10.,100.]:
      ax.plot( rfac*xcirc, rfac*ycirc, ':', color='black', linewidth=1)
   #ax.plot( Rvir*xcirc, Rvir*ycirc, ':', color='black', linewidth=1)

   #papertxt = r'$\rm{Angl\'es-Alc\'azar \, et \, al.\, 2017, \, MNRAS, \, 472, \, L109}$'
   #papertxt = r'$\rm{Angl\'es}$-${\rmAlc\'azar \, et \, al.\, 2017, \, MNRAS, \, 472, \, L109}$'
   #ax.text( 0.01, 0.005, papertxt, {'color':'orange', 'fontsize':8}, transform=ax.transAxes, horizontalalignment='left', verticalalignment='bottom' )

   fig.savefig( outfile, dpi=150 )
   plt.close()
   return 1, vmaximum
















