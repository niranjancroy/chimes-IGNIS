#>>> execfile('bhfeedback_windzoomplot.py')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec
import os as os
from shutil import copyfile
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
import daa_bhlib as bhlib
from daa_constants import *

import imp
imp.reload(daa)
imp.reload(bhlib)




basedir = '/simons/scratch/dangles/FIRE/bhfeedback/'
ICdir = 'h113_HR_sn152'
sname = 'spawn_res251_s9e1v30f05_w10dt100m1000_TshVsh'
stag = 'e1f05'

snapnum = 902

simdir = basedir + ICdir + '/' + sname

outdir = './bhfeedback/allplots'


SkipIfDone = 1


# --- define rotation matrix
bhdata = daa.load_dict_from_hdf5(simdir + '/cenBH/centralBH.hdf5')
rotmatrix = bhlib.bhrotmatrix( bhdata, Rref=1., inum=0 )
"""
Rref = 1.   # kpc
ind_Rcm = np.argmin( np.abs(bhdata['radial_bins']-Rref) )
Lgas = bhdata['cum']['L_gas'][0,ind_Rcm,:]
zaxis = np.array([0.,0.,1.])
L_dir = Lgas / np.linalg.norm(Lgas)
rotax = np.cross(L_dir,zaxis)
rotangle = np.arccos(np.dot(L_dir,zaxis))
rotmatrix_z = daa.rotation_matrix( rotax, rotangle )
#rotmatrix_tilted = daa.rotation_matrix( [1.,0.,0.], 0.5 )
#rotmatrix = np.matmul( rotmatrix_tilted, rotmatrix_z )
rotmatrix = rotmatrix_z
"""

grids = {}
for i, dx in enumerate([ 100, 10, 1 ]):
  filename = 'tzgrid_ns%04d_dx%.3f' % ( snapnum, dx )
  gridfile = simdir + '/tzgrid' + '/' + filename  + '.hdf5'
  gridDone = bhlib.make_windgrid( simdir, snapnum, gridfile, dx=dx, rotmatrix=rotmatrix, center='BH', m_wind=3000., DesNgb=64, npix=512, SkipIfDone=SkipIfDone, P=-1 )
  grd = daa.load_dict_from_hdf5( gridfile )
  grids[str(i)] = grd['edgeon']
  

lw_ax = 1.5

vmax = [ 8.5, 9.5, 10.5 ]
vmin = [ 5.5, 6.5, 7.5]

figname = stag + '_windzoomplot_ns%04d.pdf' % snapnum
print '\n', figname
fxsz = 10
xyrat = 2./3.
fysz = fxsz * xyrat
ds = 0.02
fig = plt.figure( figsize=(fxsz,fysz), dpi=300 )
gs = GridSpec(2,3)
gs.update(left=ds, right=1-ds, bottom=ds, top=1-ds, hspace=0, wspace=0)
ax1 = plt.subplot(gs[:2,:2])
ax2 = plt.subplot(gs[0,2])
ax3 = plt.subplot(gs[1,2])

for i, ax in enumerate([ax1,ax2,ax3]):
   grd = grids[str(i)]
   dx = grd['dx']
   xrg = [ -dx, dx ]
   yrg = xrg
   extent = [xrg[0],xrg[1],xrg[0],xrg[1]]  
   ax.set_xlim(xrg); ax.set_ylim(xrg); ax.set_aspect('equal')
   for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
      item.set_visible(False)
   for axis in ['top','bottom','left','right']:
      ax.spines[axis].set_linewidth(lw_ax)

   cmap_gas = cm.Greys
   cmap_wind = cm.jet
   alpha_wind = 0.25

   # --- background gas density
   gas_map = grd['gas']['m'][:,:]
   #vmax = np.log10( gas_map.max() ) - 0.5
   #vmin = vmax - 3
   im_gas = ax.imshow( np.log10( gas_map ), vmin=vmin[i], vmax=vmax[i], cmap=cmap_gas, interpolation='bicubic', origin='low', extent=extent )
   # --- AGN wind
   vwmin = 3.5 
   vwmax = 6.5 
   wind_grid = grd['wind']['m'][:,:]
   wind_map = np.ma.masked_array( wind_grid, wind_grid <= wind_grid.min() )
   im_wind = ax.imshow( np.log10( wind_map ), vmin=vwmin, vmax=vwmax, cmap=cmap_wind, interpolation='bicubic', origin='low', extent=extent, alpha=alpha_wind )
   # --- dense ISM
   #dense_map = np.ma.masked_array( gas_map, np.log10(gas_map) <= 8.5 )
   #im_dense = ax.imshow( np.log10( dense_map ), vmin=vmin, vmax=vmax, cmap=cmap_gas, interpolation='bicubic', origin='low', extent=extent, alpha=0.15 )

   if i in [0,1]:
     dx2 = grids[str(i+1)]['dx']
     xrg2 = [ -dx2, dx2 ]
     yrg2 = xrg2
     x_sq = np.array( [ xrg2[0], xrg2[1], xrg2[1], xrg2[0], xrg2[0] ] )
     y_sq = np.array( [ yrg2[0], yrg2[0], yrg2[1], yrg2[1], yrg2[0] ] )
     ax.plot( x_sq, y_sq, 'k-', lw=1 )
     if i == 0:
       x_line = np.array( [ xrg2[1], xrg[1] ] )
       y_line = np.array( [ yrg2[1], yrg[1] ] )
       ax.plot( x_line, y_line, 'k-', lw=0.5 )
       x_line = np.array( [ xrg2[1], xrg[1] ] )
       y_line = np.array( [ yrg2[0], (yrg[1]+yrg[0])/2. ] )
       ax.plot( x_line, y_line, 'k-', lw=0.5 )
     else:
       x_line = np.array( [ xrg2[0], xrg[0] ] )
       y_line = np.array( [ yrg2[0], yrg[0] ] )
       ax.plot( x_line, y_line, 'k-', lw=0.5 )
       x_line = np.array( [ xrg2[1], xrg[1] ] )
       y_line = np.array( [ yrg2[0], yrg[0] ] )
       ax.plot( x_line, y_line, 'k-', lw=0.5 )

   bar_size = dx / 2.
   xleft = -bar_size/2.
   ybar = yrg[0] + 0.06*dx
   if i == 0:
     xleft += 0.5*dx
     ybar = yrg[0] + 0.03*dx
   colbar = 'black'
   ax.plot( [xleft,xleft+bar_size], [ybar,ybar], '-', color=colbar, linewidth=2 )
   if bar_size > 1:
     bartxt = r'%d kpc' % np.round(bar_size).astype(int)
   else:
     bartxt = r'%d pc' % np.round(1e3*bar_size).astype(int)
   ax.text( xleft+0.5*bar_size, ybar+0.01*dx, bartxt, {'color':colbar, 'fontsize':10}, 
            horizontalalignment='center', verticalalignment='bottom' )  

   if i == 0:
     ax.text( 0.02, 0.98, r'${\rm log}_{10}(\Sigma_{\rm gas|wind}/{\rm M}_{\odot}{\rm kpc}^{-2})$', {'color':'black', 'fontsize':14},
              transform=ax.transAxes, horizontalalignment='left', verticalalignment='top' )

   vmin_int = np.ceil(vmin[i]).astype(int)
   vmax_int = np.floor(vmax[i]).astype(int)
   ticks = np.arange(vmin_int,vmax_int+1,1)
   if i == 0:
     cax_points = [2*ds,2*ds,0.02,0.10]
   elif i == 1:
     cax_points = [0.67,0.87,0.015,0.10]
   else:
     cax_points = [0.67,0.39,0.015,0.10]
   cbar = plt.colorbar( im_gas, cax=fig.add_axes(cax_points), ticks=ticks)
   cbar.ax.tick_params(labelsize=10)
   cbar.ax.set_yticklabels(ticks.astype(str))
   if i == 0:
     cbar.ax.set_xlabel('Gas');  cbar.ax.xaxis.set_label_position('top');  cbar.ax.xaxis.label.set_fontsize(10)

vmin_int = np.ceil(vwmin).astype(int)
vmax_int = np.floor(vwmax).astype(int)
ticks = np.arange(vmin_int,vmax_int+1,1)
cax_points = [5*ds,2*ds,0.02,0.10]
cbar = plt.colorbar( im_wind, cax=fig.add_axes(cax_points), ticks=ticks)
cbar.ax.tick_params(labelsize=10)
cbar.ax.set_yticklabels(ticks.astype(str))
cbar.ax.set_xlabel('Wind');  cbar.ax.xaxis.set_label_position('top');  cbar.ax.xaxis.label.set_fontsize(10)



fig.savefig( outdir + '/' + figname )
plt.close()






print '\n\nDone!\n'













