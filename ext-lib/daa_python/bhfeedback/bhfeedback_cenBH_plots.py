# run as >>> execfile('bhfeedback_cenBH_plots.py')

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

import daa_lib as daa
from daa_constants import *

import imp
imp.reload(daa)



#######################################
 ### DEFINE SIMULATION DIRECTORIES ###
#######################################

IC_list = [ #'sn155',
            #'sn214',
            'sn277'
          ]

sim_list = [ 'nof_s8e1_n128',
             #'spawn_s7e1v30f05_w100dt1000m100_TshVsh',
             'spawn_s8e1v30f05_w100dt1000m100_TshVsh',
             'spawn_s9e1v30f05_w100dt1000m100_TshVsh',
             'spawn_s9e1v30f05_w100dt1000m100_TshVsh_MERGE',
             'spawn_s9e1v30f05_w10dt100m1000_TshVsh',
             'spawn_s9e1v30f05_w100dt1000m100',
             'spawn_s9e1v30f033_w100dt1000m100_TshVsh',
             'spawn_s9e1v30f009_w100dt1000m100_TshVsh',
             'kick_s9e1v15000p2_n128',
             'kick_s9e1v1500p20_n128',
             'bal_s9e1v30f05_n128',
             #'TEST_SnapOrig_nof_s8e1_n128',
             #'TEST_GizmoOrig_nof_s8e1_n128',
             #'TEST_noMetDiff_nof_s8e1_n128',
             'kick_COL_s9e1v15000p2_n128',
             'kick_COL_s9e1v1500p20_n128'
           ]

sim_col = ['black','blueviolet','darkorange', 'forestgreen', 'gray', 'red', 'yellow', 'greenyellow', 'brown', 'chocolate', 'peru' ]
#sim_lw = [1,1,1]

basedir = '/simons/scratch/dangles/FIRE/bhfeedback/h113_HR_'




skipreading = 0



outdir = './bhfeedback/cenBHplot'
if not os.path.exists(outdir):
  try:
    os.mkdir(outdir)
  except OSError:
    try:
      os.makedirs(outdir)
    except OSError:
      print '...could not create output directory: ' + outdir


if skipreading==0:

 sim = {}
 for ICdir in IC_list:
    sim[ICdir] = {}
    for sname in sim_list:
       bh_file = basedir + ICdir + '/' + sname + '/cenBH/centralBH.hdf5'   
       print '...reading ' + bh_file
       bh = daa.load_dict_from_hdf5( bh_file )
       sim[ICdir][sname] = bh

       # --- add star formation history
       #sfh_file = glob( basedir + ICdir + '/' + sname + '/cenBH/SFH_*.hdf5' )[0]
       #sfh = daa.load_dict_from_hdf5( sfh_file )
       #sim[ICdir][sname]['sfh'] = sfh

 # --- read galaxy files
 gal_file = '/simons/scratch/dangles/FIRE/bhfeedback/h113_HR_galdata_fbh_00000.hdf5'
 print '...reading ' + gal_file
 gal = daa.load_dict_from_hdf5( gal_file )
 SFR_orig = gal['0.1Rvir']['SFR']
 t_orig = 1e3 * util.age_of_universe( gal['redshift'], h=0.697, Omega_M=0.2821 )


# ---- assuming same radial bins for ALL runs
radial_bins = sim[IC_list[0]][sim_list[0]]['radial_bins'][:]
x = np.log10( radial_bins )
x = ( x[:-1] + x[1:] ) / 2.
area = np.pi * (radial_bins*1e3)**2   # in pc^2
dA = area[1:] - area[:-1]
volume = (3./4.) * np.pi * (radial_bins*1e3)**3   # in pc^3
dV = volume[1:] - volume[:-1]

itime = np.array([ 0, 9, 19, 29, 39, 49, 59, 69, 79, 89, 99, 199, 299, 399, 499, 599, 699, 799, 899, 999, 1049, 1099, 1149, 1199, 1249, 1299, 1349, 1399 ])
#itime = np.array([ 0, 9, 59, 99, 499, 999, 1049, 1099, 1149, 1199, 1249, 1299, 1349, 1399 ])



# --- plots ---

left = 0.15
right = 0.95
bottom = 0.14
top = 0.95
lw = 2

rasterized = False


vmin = np.log10(1.)
vmax = np.log10(51.)
t_ticks = np.array([0,2,10,50])
ticks = np.log10(1.+t_ticks)
ticksL = t_ticks.astype(str)
cNorm  = cl.Normalize(vmin=vmin, vmax=vmax)
sMap = plt.cm.ScalarMappable(norm=cNorm, cmap=plt.cm.get_cmap('Spectral'))



for ICdir in IC_list:
 for stag in sim_list:
 
  thesim = sim[ICdir][stag]

  ylim = [-4,4]
  xlim = [-2,1]
  figname = 'SigmaGasV_vs_r_' + ICdir + '_' + stag + '.pdf'
  print '\n'+figname
  fig = plt.figure( figsize=(6,5), dpi=600 )
  ax = fig.add_subplot(111)
  fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
  for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
     item.set_fontsize(18)
  ax.set_ylabel(r'${\rm log}_{10}(\Sigma/{\rm M}_{\odot}{\rm pc}^{-3})$');  ax.set_ylim(ylim)
  ax.set_xlabel(r'${\rm log}_{10}(R \/ / \/ {\rm kpc})$');  ax.set_xlim(xlim) 
  t = thesim['time'][:] - thesim['time'][0]
  nsnap = t.size
  sMap.set_array(t)
  #for i in np.arange(nsnap):
  for i in itime[itime<nsnap]:
     Mgas = thesim['sph']['GasMass'][i,:]
     ind = np.where( thesim['sph']['Ngas'][i,:] > 0 )[0]
     xgas = x[ind]
     sigma_gas = np.log10(Mgas[ind]/dV[ind])   
     ax.plot( xgas, sigma_gas, '-', color=sMap.to_rgba(np.log10(1+t[i])), lw=2 )
  cbar = fig.colorbar( sMap, cax=fig.add_axes([0.78,0.7,0.05,0.2]), label=r'$t ({\rm Myr})$', ticks=ticks )
  cbar.ax.set_yticklabels(ticksL)
  fig.savefig( outdir + '/' + figname )
  plt.close()

  """
  ylim = [-4,4]
  xlim = [-2,1]
  figname = 'SigmaGas_vs_r_' + ICdir + '_' + stag + '.pdf'
  print '\n'+figname
  fig = plt.figure( figsize=(6,5), dpi=600 )
  ax = fig.add_subplot(111)
  fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
  for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
     item.set_fontsize(18)
  ax.set_ylabel(r'${\rm log}_{10}(\Sigma/{\rm M}_{\odot}{\rm pc}^{-2})$');  #ax.set_ylim(ylim)
  ax.set_xlabel(r'${\rm log}_{10}(R \/ / \/ {\rm kpc})$');  ax.set_xlim(xlim)
  t = thesim['time'][:] - thesim['time'][0]
  nsnap = t.size
  sMap.set_array(t)
  #for i in np.arange(nsnap):
  for i in itime[itime<nsnap]:
     Mgas = thesim['cyl']['1kpc']['GasMass'][i,:]
     ind = np.where( thesim['cyl']['1kpc']['Ngas'][i,:] > 0 )[0]
     xgas = x[ind]
     sigma_gas = np.log10(Mgas[ind]/dA[ind])
     ax.plot( xgas, sigma_gas, '-', color=sMap.to_rgba(np.log10(1+t[i])), lw=2 )
  cbar = fig.colorbar( sMap, cax=fig.add_axes([0.78,0.7,0.05,0.2]), label=r'$t ({\rm Myr})$', ticks=ticks )
  cbar.ax.set_yticklabels(ticksL)
  fig.savefig( outdir + '/' + figname )
  plt.close()
  """


  ylim = [-4,4]
  xlim = [-2,1]
  figname = 'SigmaSFR_vs_r_' + ICdir + '_' + stag + '.pdf'
  print '\n'+figname
  fig = plt.figure( figsize=(6,5), dpi=600 )
  ax = fig.add_subplot(111)
  fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
  for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
     item.set_fontsize(18)
  ax.set_ylabel(r'${\rm log}_{10}(\Sigma_{\rm SFR}/{\rm M}_{\odot}{\rm yr}^{-1}{\rm pc}^{-3})$');  #ax.set_ylim(ylim)
  ax.set_xlabel(r'${\rm log}_{10}(R \/ / \/ {\rm kpc})$');  ax.set_xlim(xlim)
  t = thesim['time'][:] - thesim['time'][0]
  nsnap = t.size
  sMap.set_array(t)
  #for i in np.arange(nsnap):
  for i in itime[itime<nsnap]:
     SFR = thesim['sph']['SFR'][i,:]
     ind = np.where( thesim['sph']['Ngas'][i,:] > 0 )[0]
     xgas = x[ind]
     sigma_SFR = np.log10(SFR[ind]/dV[ind])
     ax.plot( xgas, sigma_SFR, '-', color=sMap.to_rgba(np.log10(1+t[i])), lw=2 )
  cbar = fig.colorbar( sMap, cax=fig.add_axes([0.78,0.7,0.05,0.2]), label=r'$t ({\rm Myr})$', ticks=ticks )
  cbar.ax.set_yticklabels(ticksL)
  fig.savefig( outdir + '/' + figname )
  plt.close()




ylim = [-0.5,4]
xlim = [0,np.log10(51.)]

for ICdir in IC_list:

 figname = 'SFR_vs_t_' + ICdir + '.pdf'
 print '\n'+figname
 fig = plt.figure( figsize=(6,5), dpi=600 )
 ax = fig.add_subplot(111)
 fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
 for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(18)
 ax.set_ylabel(r'${\rm log}_{10}(SFR/{\rm M}_{\odot}{\rm yr}^{-1})$');  #ax.set_ylim(ylim)
 ax.set_xlabel(r'${\rm log}_{10}(1+t/{\rm Myr})$');  #ax.set_xlim(xlim)
 for stag, scolor in zip( sim_list, sim_col ):
    thesim = sim[ICdir][stag]
    t = thesim['time'][:] - thesim['time'][0]
    #t = np.log10( thesim['time'][:] - thesim['time'][0] + 1 )
    indR = np.argmin( np.abs(radial_bins[1:] - 100.) )
    SFR = np.log10( thesim['cum']['SFR'][:,indR] )
    ax.plot( t, SFR, '-', color=scolor, lw=1, label=stag )

    # --- from stellar ages:
    #t_ages = np.log10( thesim['sfh']['time_recent'] + 1 )
    #SFR_ages = thesim['sfh']['SFR_recent']
    #SFR_smooth = util.smooth( SFR_ages/0.85, window_len=50, window='flat' )
    #ax.plot( t_ages, np.log10( SFR_smooth ), '--', color=scolor, lw=1 )

 # --- original galaxy
 dt_orig = t_orig - thesim['time'][0]
 ind_t = np.where( (dt_orig>-30.) & (dt_orig<50.) )[0]
 ax.plot( dt_orig[ind_t], np.log10( SFR_orig[ind_t] ), '-', color='red', lw=1)
 ax.plot( dt_orig[ind_t], np.log10( SFR_orig[ind_t] ), '+', color='red', markersize=10, markeredgewidth=4 )

 l = plt.legend( loc='best', fontsize=14, frameon=False, handletextpad=0.5 )
 for line, text in zip( l.get_lines(), l.get_texts() ):
    text.set_color(line.get_color())
 fig.savefig( outdir + '/' + figname )
 plt.close()





"""

figname = 'SFR_vs_t_balVSkick.pdf'
print '\n'+figname
fig = plt.figure( figsize=(6,5), dpi=600 )
ax = fig.add_subplot(111)
fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
lw = 2
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
   item.set_fontsize(18)
ax.set_ylabel(r'${\rm log}_{10}(SFR/{\rm M}_{\odot}{\rm yr}^{-1})$');  ax.set_ylim(ylim)
ax.set_xlabel(r'${\rm log}_{10}(1+t/{\rm Myr})$');  ax.set_xlim(xlim)
for stag, scolor, slw in zip( ['nof','bal_s8','kick_s8p2','kick_s8p20'], ['lightgray','limegreen','blue','red'], [6,2,2,2] ):
 t = np.log10( sim[stag]['time'][:] - sim[stag]['time'][0] + 1 )
 indR = np.where( radial_bins <= 10. )[0]
 SFR = np.log10( sim[stag]['prf']['SFR'][:,indR].sum(axis=1) )
 ax.plot( t, SFR, '-', color=scolor, lw=slw, label=stag )
l = plt.legend( loc='lower left', fontsize=14, frameon=False, handletextpad=0.5 )
for line, text in zip( l.get_lines(), l.get_texts() ):
   text.set_color(line.get_color())
fig.savefig( outdir + '/' + figname )
plt.close()




ylim = [0.01,1.02]
xlim = [0,2.2]

figname = 'Mout_vs_t_bal.pdf'
print '\n'+figname
fig = plt.figure( figsize=(6,5), dpi=600 )
ax = fig.add_subplot(111)
fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
lw = 2
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
   item.set_fontsize(18)
ax.set_ylabel(r'$M_{\rm out}/M_{\rm gas}$');  ax.set_ylim(ylim)
ax.set_xlabel(r'${\rm log}_{10}(1+t/{\rm Myr})$');  ax.set_xlim(xlim)
for stag, scolor, slw in zip( ['nof','bal_s7','bal_s8','bal_s9'], ['lightgray','blueviolet','limegreen','darkorange'], [6,2,2,2] ):
 t = np.log10( sim[stag]['time'][:] - sim[stag]['time'][0] + 1 )
 indR = np.where( radial_bins <= 10. )[0]
 Mout = sim[stag]['prf']['GasMassOutflow'][:,indR].sum(axis=1)
 Mgas = sim[stag]['prf']['GasMass'][:,indR].sum(axis=1)
 ind = np.where( sim[stag]['prf']['Ngas'][:,indR].sum(axis=1) > 5 )[0]
 ax.plot( t[ind], Mout[ind]/Mgas[ind], '-', color=scolor, lw=slw, label=stag )
l = plt.legend( loc='lower left', fontsize=14, frameon=False, handletextpad=0.5 )
for line, text in zip( l.get_lines(), l.get_texts() ):
   text.set_color(line.get_color())
fig.savefig( outdir + '/' + figname )
plt.close()

figname = 'Mout_vs_t_balVSkick.pdf'
print '\n'+figname
fig = plt.figure( figsize=(6,5), dpi=600 )
ax = fig.add_subplot(111)
fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
lw = 2
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
   item.set_fontsize(18)
ax.set_ylabel(r'$M_{\rm out}/M_{\rm gas}$');  ax.set_ylim(ylim)
ax.set_xlabel(r'${\rm log}_{10}(1+t/{\rm Myr})$');  ax.set_xlim(xlim)
for stag, scolor, slw in zip( ['nof','bal_s8','kick_s8p2','kick_s8p20'], ['lightgray','limegreen','blue','red'], [6,2,2,2] ):
 t = np.log10( sim[stag]['time'][:] - sim[stag]['time'][0] + 1 )
 indR = np.where( radial_bins <= 10. )[0]
 Mout = sim[stag]['prf']['GasMassOutflow'][:,indR].sum(axis=1)
 Mgas = sim[stag]['prf']['GasMass'][:,indR].sum(axis=1)
 ax.plot( t, Mout/Mgas, '-', color=scolor, lw=slw, label=stag )
l = plt.legend( loc='lower left', fontsize=14, frameon=False, handletextpad=0.5 )
for line, text in zip( l.get_lines(), l.get_texts() ):
   text.set_color(line.get_color())
fig.savefig( outdir + '/' + figname )
plt.close()





ylim = [1.,4.]
xlim = [0,2.2]

figname = 'Vout_vs_t_bal.pdf'
print '\n'+figname
fig = plt.figure( figsize=(6,5), dpi=600 )
ax = fig.add_subplot(111)
fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
lw = 2
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
   item.set_fontsize(18)
ax.set_ylabel(r'${\rm log}_{10}(V_{\rm out}/{\rm km\/s}^{-1})$');  ax.set_ylim(ylim)
ax.set_xlabel(r'${\rm log}_{10}(1+t/{\rm Myr})$');  ax.set_xlim(xlim)
for stag, scolor, slw in zip( ['nof','bal_s7','bal_s8','bal_s9'], ['lightgray','blueviolet','green','darkorange'], [6,2,2,2] ):
 t = np.log10( sim[stag]['time'][:] - sim[stag]['time'][0] + 1 )
 indR = np.where( radial_bins <= 10. )[0]
 Vout = sim[stag]['prf']['GasOutflowVel'][:,indR].max(axis=1)
 #Vout = sim[stag]['prf']['GasOutflowVel_perc'][:,indR,2].max(axis=1)
 ax.plot( t, np.log10(Vout), '-', color=scolor, lw=slw, label=stag )
l = plt.legend( loc='lower left', fontsize=14, frameon=False, handletextpad=0.5 )
for line, text in zip( l.get_lines(), l.get_texts() ):
   text.set_color(line.get_color())
fig.savefig( outdir + '/' + figname )
plt.close()




ylim = [10,15]
xlim = [0,2.2]

figname = 'RadialMomentum_vs_t_bal.pdf'
print '\n'+figname
fig = plt.figure( figsize=(6,5), dpi=600 )
ax = fig.add_subplot(111)
fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
lw = 2
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
   item.set_fontsize(18)
ax.set_ylabel(r'${\rm log}_{10}(P_{\rm out}/{\rm M}_{\odot}{\rm km\/s}^{-1})$');  ax.set_ylim(ylim)
ax.set_xlabel(r'${\rm log}_{10}(1+t/{\rm Myr})$');  ax.set_xlim(xlim)
for stag, scolor, slw in zip( ['nof','bal_s7','bal_s8','bal_s9'], ['lightgray','blueviolet','green','darkorange'], [6,2,2,2] ):
 t = np.log10( sim[stag]['time'][:] - sim[stag]['time'][0] + 1 )
 indR = np.where( radial_bins <= 10. )[0]
 Pout = ( sim[stag]['prf']['GasMassOutflow'][:,indR] * sim[stag]['prf']['GasOutflowVel'][:,indR] ).sum(axis=1)
 ax.plot( t, np.log10(Pout), '-', color=scolor, lw=slw, label=stag )
l = plt.legend( loc='lower left', fontsize=14, frameon=False, handletextpad=0.5 )
for line, text in zip( l.get_lines(), l.get_texts() ):
   text.set_color(line.get_color())
fig.savefig( outdir + '/' + figname )
plt.close()

"""











print '\nDone!!\n'




