
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
import gadget as g
import utilities as util
import visualization as vis
import daa_lib as daa
from daa_constants import *

import imp
imp.reload(daa)
imp.reload(util)




def make_map( simdir, snapnum, outfile, dx=0.1, L=np.array([0.,0.,1.]), tilt=-1., hsml_fac=-1., center='COM', npix=512, SkipIfDone=0, P=-1 ):

   if SkipIfDone and os.path.exists( outfile ):
    print "\nfile exists:  %s\n" % outfile
    return 1

   # --- black hole data
   bhdata = daa.load_dict_from_hdf5(simdir + '/centralBH.hdf5')
   ind_snap = np.where( bhdata['snapnum'] == snapnum)[0];  
   if ind_snap.size == 0:
     print simdir,  'snapnum not available ', snapnum, bhdata['snapnum'].max()
     return 0
   ind_snap = ind_snap[0]
   if center == 'COM':
     ind_Rcm = np.argmin( np.abs(bhdata['cum']['radius']-0.1) )    # 100 pc
     cen_pos = bhdata['cum']['COMstar_pos'][ind_snap,ind_Rcm,:] + bhdata['bh']['pos'][ind_snap,:]    # COM of the stars
   else:
     cen_pos = bhdata['bh']['pos'][ind_snap,:]    # BH position

   # --- rotation matrix ---
   zaxis = np.array([0.,0.,1.])
   L_dir = L / np.linalg.norm(L)
   rotax = np.cross(L_dir,zaxis)
   rotangle = np.arccos(np.dot(L_dir,zaxis))
   rotmatrix = daa.rotation_matrix( rotax, rotangle )
   if tilt > 0:
     rotmatrix_tilted = daa.rotation_matrix( [1.,0.,0.], tilt )
     rotmatrix = np.matmul( rotmatrix_tilted, rotmatrix )

   # --- gas ---
   if type(P) == 'dict':
     print ' ...making grid with given P:', simdir, 'snapnum, dx =', snapnum, dx
   else:
     print ' ...making grid:', simdir, '  [ snapnum, dx =', snapnum, dx, ']'
     P = g.readsnap( simdir, snapnum, 0, cosmological=1 )
   r_gas = P['p'][:,:] - cen_pos[:]
   r_gas_rot = np.tensordot( rotmatrix, r_gas[:,:], axes=(1,1) ).T
   x = r_gas_rot[:,0]
   y = r_gas_rot[:,1]
   z = r_gas_rot[:,2]
   m = P['m'][:] * UnitMass_in_Msun
   T = g.gas_temperature(P['u'],P['ne'])
   hsml = P['h'][:]

   # --- stars
   #P = g.readsnap( simdir, snapnum, 4, cosmological=1 )
   #x = P['p'][:,0] - cen_pos[0]
   #y = P['p'][:,1] - cen_pos[1]
   #z = P['p'][:,2] - cen_pos[2]
   #m = P['m'][:] * UnitMass_in_Msun
   #hsml = daa.get_particle_hsml( x, y, z, DesNgb=32 )


   grid = { 'mgas':{ 'dx':dx, 'Npart':0, 'g':np.zeros([npix,npix]) },
            'T':{ 'dx':dx, 'Npart':0, 'g':np.zeros([npix,npix]) },
            'mgas_edgeon':{ 'dx':dx, 'Npart':0, 'g':np.zeros([npix,npix]) },
            'T_edgeon':{ 'dx':dx, 'Npart':0, 'g':np.zeros([npix,npix]) },
            'bhpos':np.zeros(3)
          }

   grid['snapnum'] = snapnum
   grid['redshift'] = P['redshift']
   bhpos = bhdata['bh']['pos'][ind_snap,:] - cen_pos
   bhpos = np.tensordot( rotmatrix, bhpos, axes=(1,0) )
   grid['bhpos'][:] = bhpos

   bfr = 3.
   xrg = [ -dx, dx ]
   ind = np.where( (np.abs(x) <= xrg[1]*bfr) & (np.abs(y) <= xrg[1]*bfr) & (np.abs(z) <= xrg[1]*bfr) )[0]
  
   if ind.size <= 1:
     print '--> no particles found in grid size <--'
     return 0
   x_grid = x[ind] 
   y_grid = y[ind] 
   z_grid = z[ind]
   m_grid = m[ind]
   T_grid = T[ind]
   hsml_grid = hsml[ind]

   if hsml_fac > 0.:
     hsml_grid *= hsml_fac
   if dx < 0.09:
     hsml_grid *= 1.5

   grid['mgas']['Npart'] = ind.size
   grid['T']['Npart'] = ind.size
   m_sMap, Tm_sMap = daa.make_2Dgrid( x_grid, y_grid, xrg, weight1=m_grid, weight2=T_grid*m_grid, hsml=hsml_grid, pixels=npix )
   Mgas_map, Tgas_map = daa.clip_2Dgrid( m_sMap, map2=Tm_sMap )
   grid['mgas']['g'][:,:] = Mgas_map.T
   grid['T']['g'][:,:] = Tgas_map.T
  
   grid['mgas_edgeon']['Npart'] = ind.size
   grid['T_edgeon']['Npart'] = ind.size
   m_sMap, Tm_sMap = daa.make_2Dgrid( x_grid, z_grid, xrg, weight1=m_grid, weight2=T_grid*m_grid, hsml=hsml_grid, pixels=npix )
   Mgas_map, Tgas_map = daa.clip_2Dgrid( m_sMap, map2=Tm_sMap )
   grid['mgas_edgeon']['g'][:,:] = Mgas_map.T
   grid['T_edgeon']['g'][:,:] = Tgas_map.T   


   ############################
    ### WRITE DATA TO FILE ###
   ############################

   if os.path.isfile(outfile):
     os.remove(outfile)
   file = h5py.File(outfile, 'w')
   for grname in grid.keys():
      if type(grid[grname]) is dict:
        group = file.create_group( grname )
        for keyname in grid[grname].keys():
           group.create_dataset(keyname, data=grid[grname][keyname])
      else:
        file.create_dataset(grname, data=grid[grname])
   file.close()
   return 1




def plot_map( gridfile, outfile, varstr='mgas', vmax=-999, galdata='', bhpos=0, rotmatrix=0, small=1, text=1, dxsqr=0, cbar=1 ):
   try:
     grd = daa.load_dict_from_hdf5( gridfile )
     print '\n', gridfile, ' --> ', outfile, '\n'
   except IOError:
     print 'grid not found: ', gridfile
     return 0, vmax
   fig = plt.figure( figsize=(5.,5.), dpi=300 )
   fig.subplots_adjust(left=0., right=1., bottom=0., top=1., hspace=0., wspace=0.)
   ax = fig.add_subplot(111)

   if 'faceon' in grd.keys():
     grd = grd['faceon']

   vmaximum = np.log10( grd[varstr]['g'][:,:].max() ) - 0.5
   if vmax == -999:
     vmax = vmaximum
   vmin = vmax - 3
   dx = grd[varstr]['dx']
   xrg = [ -dx, dx ]
   extent = [xrg[0],xrg[1],xrg[0],xrg[1]]
   ax.set_xlim(xrg); ax.set_ylim(xrg);
   ax.set_aspect('equal'); ax.axis('off')
   for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
      item.set_visible(False)
   cmap = cm.viridis
   #cmap = cm.Spectral
   #cmap = cm.afmhot
   img = ax.imshow( np.log10( grd[varstr]['g'][:,:] ), vmin=vmin, vmax=vmax, cmap=cmap, interpolation='bicubic', origin='low', extent=extent )

   if len(galdata) > 0:
     acirc = np.linspace(0,2*np.pi,100)
     xcirc = np.cos(acirc)
     ycirc = np.sin(acirc)
     R0 = 0.01
     ax.plot( R0*xcirc, R0*ycirc, '-', color='black', linewidth=1.)
     galdata = daa.load_dict_from_hdf5( galdata )
     ind_snap = np.where( galdata['snapnum'] == 156)[0]
     Rvir = galdata['halo']['Rvir'][ind_snap[0]]
     print '   --> Rvir =', Rvir, ' [kpc]'
     for rfac in [1.,0.1,0.01]:
        ax.plot( rfac*Rvir*xcirc, rfac*Rvir*ycirc, ':', color='white', linewidth=1)

   if bhpos == 1:
     lw = 0.5
     s = 20
     if util.checklen(rotmatrix) == 1:
       bhpos = grd['bhpos'][:]
     else:
       bhpos = np.tensordot( rotmatrix, grd['bhpos'][:], axes=(1,0) )
     if varstr == 'mgas_edgeon':
       ax.scatter( bhpos[0], bhpos[2], s=s, c='black', marker='+', linewidth=lw )
     else:
       ax.scatter( bhpos[0], bhpos[1], s=s, c='black', marker='+', linewidth=lw )

   if dxsqr > 0:
     x_sq = dxsqr * np.array( [ -1, 1, 1, -1, -1 ] )
     y_sq = dxsqr * np.array( [ -1, -1, 1, 1, -1 ] )
     ax.plot( x_sq, y_sq, '-', color='white', lw=1 )

   if cbar == 1:
     vmin_int = np.ceil(vmin).astype(int)
     vmax_int = np.floor(vmax).astype(int)
     ticks = np.arange(vmin_int,vmax_int+1,1)
     if small:
       cbar = plt.colorbar( img, cax=fig.add_axes([1-0.07,0.+0.01,0.02,0.10]), ticks=ticks)#,
       cbar.ax.tick_params(labelsize=8)
     else:
       cbar = plt.colorbar( img, cax=fig.add_axes([1-0.10,0.+0.02,0.03,0.15]), ticks=ticks)#,
       cbar.ax.tick_params(labelsize=11)
     cbar.ax.set_yticklabels(ticks.astype(str))
     for this in ['yticklabels']:
        cbytick_obj = plt.getp(cbar.ax.axes, this)
        plt.setp(cbytick_obj, color='white')
     cbar.ax.yaxis.set_tick_params(color='w')
     cbar.outline.set_edgecolor('white')

   if text:
     if small:
       fontsize = 10
     else:
       fontsize = 15
     #ax.text( 0.99, 0.99, r'${\rm log}(\Sigma_{\rm gas}/{\rm M}_{\odot}{\rm kpc}^{-2})$', {'color':'black', 'fontsize':fontsize}, bbox={'color':'white', 'pad':0.5},
     #         transform=ax.transAxes, horizontalalignment='right', verticalalignment='top' )
     ax.text( 0.01, 0.99, r'${\rm log}_{10}(\Sigma_{\rm gas}/{\rm M}_{\odot}{\rm kpc}^{-2})$', {'color':'black', 'fontsize':fontsize}, bbox={'color':'white', 'pad':0.5},
              transform=ax.transAxes, horizontalalignment='left', verticalalignment='top' )

   fig.savefig( outfile )
   plt.close()
   return 1, vmaximum







