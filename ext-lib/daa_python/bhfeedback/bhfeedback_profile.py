# run as >>> execfile('bhfeedback_profile.py')
# python bhfeedback_profile.py > log/prfBH.0.log 2>&1 &

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

basedir = '/simons/scratch/dangles/FIRE/bhfeedback/'
ICdir = 'h113_HR_sn214/'

sim_list = [ 'nof_s8e1_n256',
             'bal_s9e1v30f05_n256',
             'bal_s8e1v30f05_n256',
             'bal_s7e1v30f05_n256',
             'kick_s8e1v15000p2_n256',
             'kick_s8e1v1500p20_n256',
             'kick_s7e1v15000p2_n256',
             'kick_s7e1v1500p20_n256' ]

sim_tag = [ 'nof',
            'bal_s9',
            'bal_s8',
            'bal_s7',
            'kick_s8p2',
            'kick_s8p20',
            'kick_s7p2',
            'kick_s7p20' ]

nsim = len(sim_list)


skipreading = 1



outdir = './bhfeedback/'+ ICdir + 'prf'
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

 for sname, stag in zip(sim_list,sim_tag):

    simdir = basedir + ICdir + sname
    bh_file = simdir + '/bh_rbins.hdf5'   
 
    print '...reading ' + bh_file
    bh = daa.load_dict_from_hdf5( bh_file )
    sim[stag] = bh





# ---- assuming same radial bins for ALL runs
radial_bins = sim[sim_tag[0]]['prf']['radial_bins'][:]
x = np.log10( radial_bins )
x = ( x[:-1] + x[1:] ) / 2.
area = np.pi * (radial_bins*1e3)**2   # in pc^2
da = area[1:] - area[:-1]
volume = (3./4.) * np.pi * (radial_bins*1e3)**3   # in pc^3
dV = volume[1:] - volume[:-1]


# --- plots ---

left = 0.15
right = 0.95
bottom = 0.14
top = 0.95

rasterized = False

ylim = [-6,3]
xlim = [-2,1]

vmin = np.log10(1.)
vmax = np.log10(100.)
ticks = np.array([0.0,0.5,1.0,1.5,2.0])
cNorm  = cl.Normalize(vmin=vmin, vmax=vmax)
sMap = plt.cm.ScalarMappable(norm=cNorm, cmap=plt.cm.get_cmap('Spectral'))


figname = 'SigmaGas_vs_r_bal.pdf'
print '\n'+figname
fig = plt.figure( figsize=(6,5), dpi=600 )
ax = fig.add_subplot(111)
fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
lw = 2
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
   item.set_fontsize(18)
ax.set_ylabel(r'${\rm log}_{10}(\Sigma/{\rm M}_{\odot}{\rm pc}^{-3})$');  ax.set_ylim(ylim)
ax.set_xlabel(r'${\rm log}_{10}(R \/ / \/ {\rm kpc})$');  ax.set_xlim(xlim)
#for stag in ['bal_s9','bal_s8','bal_s7']:
for stag in ['bal_s9']:
 t = sim[stag]['time'][:] - sim[stag]['time'][0]
 nsnap = t.size
 sMap.set_array(t)
 for i in np.arange(nsnap):
   Mgas = sim[stag]['prf']['GasMass'][i,:]
   ind = np.where( sim[stag]['prf']['Ngas'][i,:] > 5 )[0]
   xgas = x[ind]
   sigma_gas = np.log10(Mgas[ind]/dV[ind])   
   #ax.plot( xgas, sigma_gas, '-', color=sMap.to_rgba(t[i]), lw=2 )
   #sigma_gas = np.log10(Mgas/dV)
   #ax.plot( x, sigma_gas, '-', color=sMap.to_rgba(t[i]), lw=2 )
   ax.plot( xgas, sigma_gas, '-', color=sMap.to_rgba(np.log10(1+t[i])), lw=2 )
cbar = fig.colorbar( sMap, cax=fig.add_axes([0.78,0.7,0.05,0.2]), label=r'${\rm log}_{10}(1+t/{\rm Myr})$', ticks=ticks )
cbar.ax.set_yticklabels(ticks.astype(str))
fig.savefig( outdir + '/' + figname )
plt.close()



"""
figname = 'SigmaGas_vs_r_bal_4p.pdf'
print '\n'+figname
fig = plt.figure( figsize=(6,5), dpi=600 )
fig.subplots_adjust(left=0.12, right=0.96, bottom=0.11, top=0.96, hspace=0, wspace=0)
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222, sharey=ax1)
ax3 = fig.add_subplot(223, sharex=ax1)
ax4 = fig.add_subplot(224, sharex=ax2, sharey=ax3)
lw = 0.5
for i, ax, stag in zip( range(4), [ax1,ax2,ax3,ax4], ['nof','bal_s7','bal_s8','bal_s9'] ):
 if rasterized:
   ax.set_rasterization_zorder(1)
 for item in ([ax.title] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(11)
 for item in [ax.xaxis.label, ax.yaxis.label]:
    item.set_fontsize(15)
 ax.set_ylim(ylim);  ax.set_yticks(np.arange(-5,3))
 ax.set_xlim(xlim);  ax.set_xticks(np.arange(-2,1,0.5))
 ax.text( 0.98, 0.98, stag, {'color':'black', 'fontsize':14}, bbox={'color':'white', 'pad':0.4},
         transform=ax.transAxes, horizontalalignment='right', verticalalignment='top' )
 if i in [0,2]:
   ax.set_ylabel(r'${\rm log}_{10}(\Sigma/{\rm M}_{\odot}{\rm pc}^{-3})$')
 else:
   for label in ax.get_yticklabels():
      label.set_visible(False)
 if i in [2,3]:
   ax.set_xlabel(r'${\rm log}_{10}(R \/ / \/ {\rm kpc})$')
 else:
   for label in ax.get_xticklabels():
      label.set_visible(False)
 t = sim[stag]['time'][:] - sim[stag]['time'][0]
 nsnap = t.size
 sMap.set_array(t);
 for i in np.arange(nsnap):
   Mgas = sim[stag]['prf']['GasMass'][i,:]
   ind = np.where( sim[stag]['prf']['Ngas'][i,:] > 5 )[0]
   xgas = x[ind]
   sigma_gas = np.log10(Mgas[ind]/dV[ind])
   ax.plot( xgas, sigma_gas, '-', color=sMap.to_rgba(np.log10(1+t[i])), lw=2 )
cbar = fig.colorbar( sMap, cax=fig.add_axes([0.14,0.57,0.05,0.2]), label=r'${\rm log}_{10}(1+t/{\rm Myr})$', ticks=ticks )
cbar.ax.set_yticklabels(ticks.astype(str))
for item in cbar.ax.get_yticklabels():
   item.set_fontsize(9)
fig.savefig( outdir + '/' + figname )
plt.close()


figname = 'SigmaGas_vs_r_balVSkick_4p.pdf'
print '\n'+figname
fig = plt.figure( figsize=(6,5), dpi=600 )
fig.subplots_adjust(left=0.12, right=0.96, bottom=0.11, top=0.96, hspace=0, wspace=0)
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222, sharey=ax1)
ax3 = fig.add_subplot(223, sharex=ax1)
ax4 = fig.add_subplot(224, sharex=ax2, sharey=ax3)
lw = 0.5
for i, ax, stag in zip( range(4), [ax1,ax2,ax3,ax4], ['nof','bal_s8','kick_s8p2','kick_s8p20'] ):
 if rasterized:
   ax.set_rasterization_zorder(1)
 for item in ([ax.title] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(11)
 for item in [ax.xaxis.label, ax.yaxis.label]:
    item.set_fontsize(15)
 ax.set_ylim(ylim);  ax.set_yticks(np.arange(-5,3))
 ax.set_xlim(xlim);  ax.set_xticks(np.arange(-2,1,0.5))
 ax.text( 0.98, 0.98, stag, {'color':'black', 'fontsize':14}, bbox={'color':'white', 'pad':0.4},
         transform=ax.transAxes, horizontalalignment='right', verticalalignment='top' )
 if i in [0,2]:
   ax.set_ylabel(r'${\rm log}_{10}(\Sigma/{\rm M}_{\odot}{\rm pc}^{-3})$')
 else:
   for label in ax.get_yticklabels():
      label.set_visible(False)
 if i in [2,3]:
   ax.set_xlabel(r'${\rm log}_{10}(R \/ / \/ {\rm kpc})$')
 else:
   for label in ax.get_xticklabels():
      label.set_visible(False)
 t = sim[stag]['time'][:] - sim[stag]['time'][0]
 nsnap = t.size
 sMap.set_array(t);
 for i in np.arange(nsnap):
   Mgas = sim[stag]['prf']['GasMass'][i,:]
   ind = np.where( sim[stag]['prf']['Ngas'][i,:] > 5 )[0]
   xgas = x[ind]
   sigma_gas = np.log10(Mgas[ind]/dV[ind])
   ax.plot( xgas, sigma_gas, '-', color=sMap.to_rgba(np.log10(1+t[i])), lw=2 )
cbar = fig.colorbar( sMap, cax=fig.add_axes([0.14,0.57,0.05,0.2]), label=r'${\rm log}_{10}(1+t/{\rm Myr})$', ticks=ticks )
cbar.ax.set_yticklabels(ticks.astype(str))
for item in cbar.ax.get_yticklabels():
   item.set_fontsize(9)
fig.savefig( outdir + '/' + figname )
plt.close()
"""




ylim = [-1.99,3]
xlim = [0,2.2]

figname = 'SFR_vs_t_bal.pdf'
print '\n'+figname
fig = plt.figure( figsize=(6,5), dpi=600 )
ax = fig.add_subplot(111)
fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
lw = 2
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
   item.set_fontsize(18)
ax.set_ylabel(r'${\rm log}_{10}(SFR/{\rm M}_{\odot}{\rm yr}^{-1})$');  ax.set_ylim(ylim)
ax.set_xlabel(r'${\rm log}_{10}(1+t/{\rm Myr})$');  ax.set_xlim(xlim)
for stag, scolor, slw in zip( ['nof','bal_s7','bal_s8','bal_s9'], ['lightgray','blueviolet','limegreen','darkorange'], [6,2,2,2] ):
 t = np.log10( sim[stag]['time'][:] - sim[stag]['time'][0] + 1 )
 indR = np.where( radial_bins <= 10. )[0]
 SFR = np.log10( sim[stag]['prf']['SFR'][:,indR].sum(axis=1) )
 ax.plot( t, SFR, '-', color=scolor, lw=slw, label=stag )
l = plt.legend( loc='lower left', fontsize=14, frameon=False, handletextpad=0.5 )
for line, text in zip( l.get_lines(), l.get_texts() ):
   text.set_color(line.get_color())
fig.savefig( outdir + '/' + figname )
plt.close()

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













print '\nDone!!\n'




