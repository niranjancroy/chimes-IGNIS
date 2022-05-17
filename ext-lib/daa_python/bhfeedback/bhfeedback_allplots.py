# run as >>> execfile('bhfeedback_allplots.py')

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

sim_list = [ 'nof_s8e1_n128',
             'spawn_res251_s9e01v30f05_w10dt100m1000_TshVsh',
             'spawn_res251_s9e1v30f05_w10dt100m1000_TshVsh',
             'spawn_res251_s9e1v30f033_w10dt100m1000_TshVsh',
             'spawn_res251_s9e1v30f02_w10dt100m1000_TshVsh',
             #'spawn_res251_s9e1v30f011_w10dt100m1000_TshVsh',
             'spawn_res251_s9e1v30f009_w10dt100m1000_TshVsh'#,
             #'TEST_mainbranch_e1f05',
             #'TEST_mainbranch_e1f009'
           ]

tag_list = [ 'nof',
             'e01f05',
             'e1f05',
             'e1f033',
             'e1f02',
             #'e1f011',
             'e1f009'#,
             #'Te1f05',
             #'Te1f009',
           ]

facc_list = [ 1.,
              0.5,
              0.5,
              0.33333333,
              0.2,
              #0.11111111,
              0.09090909#,
              #0.5,
              #0.09090909
            ]

h = 0.697
Mbh = 7e8 / h # Msun
Mdot_e1 = daa.MdotEddington( Mbh )
Mdot_e01 = 0.1 * Mdot_e1       # Msun/year

mdot_list = [ 0.,
              Mdot_e01,
              Mdot_e1,
              Mdot_e1,
              Mdot_e1,
              Mdot_e1,
              Mdot_e1#,
              #Mdot_e1,
              #Mdot_e1
            ]


sim_col = ['gray','blueviolet','dodgerblue','forestgreen','greenyellow','darkorange', 'red' ]

basedir = '/simons/scratch/dangles/FIRE/bhfeedback/h113_HR_sn152/'




skipreading = 1



outdir = './bhfeedback/allplots'
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
 for sname, stag in zip(sim_list,tag_list):
    bh_file = basedir  + sname + '/cenBH/centralBH.hdf5'   
    print '...reading ' + bh_file
    bh = daa.load_dict_from_hdf5( bh_file )
    sim[stag] = bh
 time_start = sim[tag_list[1]]['time'][0] 
 redshift_start = sim[tag_list[1]]['redshift'][0]

 # ---- assuming same radial bins for ALL runs
 radial_bins = sim[tag_list[0]]['radial_bins'][:]
 x = np.log10( radial_bins )
 x = ( x[:-1] + x[1:] ) / 2.
 nbin = x.size
 area = np.pi * (radial_bins)**2   # in kpc^2
 dA = area[1:] - area[:-1]
 volume = (3./4.) * np.pi * (radial_bins)**3   # in kpc^3
 dV = volume[1:] - volume[:-1]



# ---
fsz_axis = 16
fsz_num = 14
lw_ax = 1.5
lw = 1.5

left = 0.15
right = 0.96
bottom = 0.14
top = 0.96

tmin = 0.
tmax = 54.
t_ticks = np.array([0,20,40,60])
ticks = t_ticks
ticksL = t_ticks.astype(str)
cNorm  = cl.Normalize(vmin=tmin, vmax=tmax)
sMap = plt.cm.ScalarMappable(norm=cNorm, cmap=plt.cm.get_cmap('Spectral'))
# ---


###############################
 ### original galaxy plots ###
###############################
"""
figname = 'gal_vs_z.pdf'
print '\n'+figname
fig = plt.figure( figsize=(5,7), dpi=300 )
fig.subplots_adjust(left=0.16, right=0.98, bottom=0.08, top=0.94, hspace=0.0)
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312, sharex=ax1)
ax3 = fig.add_subplot(313, sharex=ax1)
for ax in ax1, ax2, ax3:
   ax.patch.set_visible(False)
   for item in (ax.get_xticklabels() + ax.get_yticklabels()):
      item.set_fontsize(fsz_num)
   for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
      item.set_fontsize(fsz_axis)
   for axis in ['top','bottom','left','right']:
      ax.spines[axis].set_linewidth(lw_ax)
   ax.xaxis.set_tick_params(width=lw_ax)
   ax.yaxis.set_tick_params(width=lw_ax)
   if ax != ax3:
     for label in ax.get_xticklabels():
        label.set_visible(False)
lw = 1
color = 'black'
color_gas = 'dodgerblue'
color_star = 'darkorange'
color_gal = 'limegreen' #'m' #'plum'
color_out = 'm' #'limegreen' #'green'

gal_file = '/simons/scratch/dangles/FIRE/bhfeedback/h113_HR_galdata_fbh_00000.hdf5'
print '...reading ' + gal_file
gal = daa.load_dict_from_hdf5( gal_file )
redshift_orig = gal['redshift']
t_orig = 1e3 * util.age_of_universe( redshift_orig, h=0.697, Omega_M=0.2821 )
Mstar_orig = gal['0.1Rvir']['StarMass']
SFR_orig = gal['0.1Rvir']['SFR']
rtag = '1kpc'
Mtot_orig = gal[rtag]['GasMass'] + gal[rtag]['StarMass'] + gal[rtag]['DarkMass']
ind_m = np.where( Mstar_orig > 0)[0]
Vcirc_orig = np.sqrt( G_UNIV * (Mtot_orig[ind_m]*MSUN) / (gal[rtag]['R0'][ind_m]*CM_PER_KPC) ) / CM_PER_KM     # km/s
ind_start = np.argmin( np.abs(redshift_orig - redshift_start) )
Mstar_start = Mstar_orig[ind_start]
SFR_start = SFR_orig[ind_start]
Vcirc_start = Vcirc_orig[ind_start]
print 'redshift_start =', redshift_start, '  Mstar_start =', Mstar_start, '  SFR_start =', SFR_start, '  Vcirc_start =', Vcirc_start

redshift_orig = np.log10( 1 + redshift_orig[ind_m] )
xrg = [ np.log10( 1. + 0.999 ), np.log10( 1. + 6.1 ) ]

axtw = ax1.twiny(); axtw.set_xlabel(r'$z$');  axtw.set_xlim(xrg[0],xrg[1])
ztmp = np.arange(11)
logztmp = np.log10(1+ztmp)
ztmp = ztmp[ (logztmp>xrg[0]) & (logztmp<xrg[1]) ]
axtw.set_xticks( np.log10(1+ztmp) ); axtw.set_xticklabels( ztmp.astype('str') );
for item in (axtw.get_xticklabels()):
   item.set_fontsize(fsz_num-4)
for item in ([axtw.xaxis.label]):
   item.set_fontsize(fsz_axis-2)

ax1.set_ylabel(r'${\rm log}_{10}(M_{\rm star}/{\rm M}_{\odot})$');  ax1.set_ylim(7.5,11.5);  ax1.set_yticks([8,9,10,11]);  ax1.set_xlim(xrg[0],xrg[1])
ax1.plot( redshift_orig, np.log10( Mstar_orig[ind_m] ), '-', color=color, lw=lw )

ax2.set_ylabel(r'${\rm log}_{10}(SFR/{\rm M}_{\odot}{\rm yr}^{-1})$');  ax2.set_ylim(-1.5,3.5);  ax2.set_yticks([-1,0,1,2,3])
ax2.plot( redshift_orig, np.log10( SFR_orig[ind_m] ), '-', color=color, lw=lw )

ax3.set_ylabel(r'$V_{\rm circ} ({\rm km\/s}^{-1})$');  ax3.set_ylim([0,900]);  ax3.set_yticks([200,400,600,800])
ax3.set_xlabel(r'${\rm log}_{10}( 1 + z )$')
ax3.text( 0.05, 0.08, r'$r < 1 {\rm kpc}$', {'color':'black', 'fontsize':fsz_axis}, transform=ax3.transAxes, horizontalalignment='left' )
ax3.plot( redshift_orig, Vcirc_orig , '-', color=color, lw=lw )

#z_int = np.log10(1+ np.array([sim['nof']['redshift'][0],sim['nof']['redshift'][-1]]) )
z_int = np.log10(1+ np.array([2.34,2.23]) )
z_this = np.array([redshift_start,redshift_start])
for ax in ax1, ax2, ax3:
   y_tmp = ax.get_ylim()
   ax.plot( np.log10(1+z_this), y_tmp, ':', color='black', lw=1)
   ax.fill_between(z_int, [y_tmp[0],y_tmp[0]], [y_tmp[1],y_tmp[1]], color='limegreen', zorder=10, alpha=0.5 )

fig.savefig( outdir + '/' + figname )
plt.close()



#################################
 ### single simulation plots ###
#################################

left = 0.15
right = 0.96
bottom = 0.14
top = 0.96

tmin = 0.
tmax = 54.
t_ticks = np.array([0,20,40,60])
ticks = t_ticks
ticksL = t_ticks.astype(str)
cNorm  = cl.Normalize(vmin=tmin, vmax=tmax)
sMap = plt.cm.ScalarMappable(norm=cNorm, cmap=plt.cm.get_cmap('Spectral'))

for stag in ['nof', 'e1f009']:
 
 thesim = sim[stag]

 xlim = [-2,0]
 figname = 'SigmaStar_vs_r_' + stag + '.pdf'
 print '\n'+figname
 fig = plt.figure( figsize=(5,4), dpi=600 )
 ax = fig.add_subplot(111)
 fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
 for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(fsz_num)
 for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
    item.set_fontsize(fsz_axis)
 for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(lw_ax)
 ax.xaxis.set_tick_params(width=lw_ax)
 ax.yaxis.set_tick_params(width=lw_ax)
 ylim = [9,13]
 ax.set_ylabel(r'${\rm log}_{10}(\Sigma_{\rm star}/{\rm M}_{\odot}{\rm kpc}^{-2})$');  ax.set_ylim(ylim);  ax.set_yticks([9,10,11,12,13])
 ax.set_xlabel(r'${\rm log}_{10}(R \/ / \/ {\rm kpc})$');  ax.set_xlim(xlim)
 t = thesim['time'][:] - time_start
 nsnap = t.size
 sMap.set_array(t)
 for i in np.arange(nsnap):
    if (t[i]<tmin) or (t[i]>tmax):
      continue
    Mstar = thesim['cyl']['1kpc']['StarMass'][i,:]
    ind = np.where( thesim['cyl']['1kpc']['Nstar'][i,:] > 0 )[0]
    xstar = x[ind]
    sigma_star = np.log10(Mstar[ind]/dA[ind])
    ax.plot( xstar, sigma_star, '-', color=sMap.to_rgba(t[i]), lw=lw )
 cbar = fig.colorbar( sMap, cax=fig.add_axes([0.78,0.7,0.05,0.2]), label=r'$t ({\rm Myr})$', ticks=ticks )
 cbar.ax.set_yticklabels(ticksL)
 fig.savefig( outdir + '/' + figname )
 plt.close()

 
 figname = 'SigmaGas_vs_r_' + stag + '.pdf'
 print '\n'+figname
 fig = plt.figure( figsize=(5,4), dpi=600 )
 ax = fig.add_subplot(111)
 fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
 for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(fsz_num)
 for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
    item.set_fontsize(fsz_axis)
 for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(lw_ax)
 ax.xaxis.set_tick_params(width=lw_ax)
 ax.yaxis.set_tick_params(width=lw_ax)
 ylim = [7,11]
 ax.set_ylabel(r'${\rm log}_{10}(\Sigma_{\rm gas}/{\rm M}_{\odot}{\rm kpc}^{-2})$');  ax.set_ylim(ylim);  ax.set_yticks([7,8,9,10,11])
 ax.set_xlabel(r'${\rm log}_{10}(R \/ / \/ {\rm kpc})$');  ax.set_xlim(xlim)
 t = thesim['time'][:] - time_start
 nsnap = t.size
 sMap.set_array(t)
 for i in np.arange(nsnap):
    if (t[i]<tmin) or (t[i]>tmax):
      continue
    #Mgas = thesim['cyl']['1kpc']['GasMass'][i,:]
    #Ngas_inbin = thesim['cyl']['1kpc']['Ngas'][i,:]
    Mgas = thesim['sph']['GasMass'][i,:]
    Ngas_inbin = thesim['sph']['Ngas'][i,:]
    ind = np.where( Ngas_inbin > 0 )[0]
    xgas = x[ind]
    sigma_gas = np.log10(Mgas[ind]/dA[ind])
    ax.plot( xgas, sigma_gas, '-', color=sMap.to_rgba(t[i]), lw=lw )
 cbar = fig.colorbar( sMap, cax=fig.add_axes([0.78,0.7,0.05,0.2]), label=r'$t ({\rm Myr})$', ticks=ticks )
 cbar.ax.set_yticklabels(ticksL)
 fig.savefig( outdir + '/' + figname )
 plt.close()
"""


###################################
 ### feedback comparison plots ###
##################################

nof_color = 'lightgray'
nof_lw = 3

xlim = [-10.,tmax]

figname = 'SFR_vs_t.pdf'
print '\n'+figname
fig = plt.figure( figsize=(5,4), dpi=600 )
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.17, right=right, bottom=bottom, top=top)
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
   item.set_fontsize(fsz_num)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
   item.set_fontsize(fsz_axis)
for axis in ['top','bottom','left','right']:
   ax.spines[axis].set_linewidth(lw_ax)
ax.xaxis.set_tick_params(width=lw_ax)
ax.yaxis.set_tick_params(width=lw_ax)
ax.set_ylabel(r'${\rm log}_{10}(SFR/{\rm M}_{\odot}{\rm yr}^{-1})$');  #ax.set_ylim(ylim)
ax.set_xlabel(r'${\rm Time \/ (Myr)}$');  ax.set_xlim(xlim)
for stag, scolor in zip( tag_list, sim_col ):
   thesim = sim[stag]
   t = thesim['time'][:] - time_start
   rcut = 10.
   indR = np.argmin( np.abs(radial_bins[1:] - rcut) )
   SFR = np.log10( thesim['cum']['SFR'][:,indR] )
   if stag == 'nof':
     ax.plot( t, SFR, '-', color=nof_color, lw=nof_lw, label=stag )
   else:
     ax.plot( t, SFR, '-', color=scolor, lw=lw, label=stag )
#l = plt.legend( loc='best', fontsize=14, frameon=False, handletextpad=0.5 )
#for line, text in zip( l.get_lines(), l.get_texts() ):
#   text.set_color(line.get_color())
fig.savefig( outdir + '/' + figname )
plt.close()



figname = 'Vcirc_vs_t.pdf'
print '\n'+figname
fig = plt.figure( figsize=(5,4), dpi=600 )
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.17, right=right, bottom=bottom, top=top)
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
   item.set_fontsize(fsz_num)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
   item.set_fontsize(fsz_axis)
for axis in ['top','bottom','left','right']:
   ax.spines[axis].set_linewidth(lw_ax)
ax.xaxis.set_tick_params(width=lw_ax)
ax.yaxis.set_tick_params(width=lw_ax)
ax.set_ylabel(r'$V_{\rm circ} ({\rm km\/s}^{-1})$');  #ax.set_ylim([0,900]);  ax.set_yticks([200,400,600,800])
ax.set_xlabel(r'${\rm Time \/ (Myr)}$');  ax.set_xlim(xlim)
for stag, scolor in zip( tag_list, sim_col ):
   thesim = sim[stag]
   t = thesim['time'][:] - time_start
   rcut = 1. #kpc
   indR = np.argmin( np.abs(radial_bins[1:] - rcut) )
   #Vcirc = thesim['cum']['Vcirc'][:,indR]
   Mbh = 1e9
   Mtot = thesim['cum']['DarkMass'][:,indR] + thesim['cum']['StarMass'][:,indR] + thesim['cum']['GasMass'][:,indR] + Mbh
   R0 = radial_bins[1:][indR]
   Vcirc = np.sqrt( G_UNIV * (Mtot*MSUN) / (R0*CM_PER_KPC) ) / CM_PER_KM     # km/s
   if stag == 'nof':
     ax.plot( t, Vcirc, '-', color=nof_color, lw=nof_lw, label=stag )
   else:
     ax.plot( t, Vcirc, '-', color=scolor, lw=lw, label=stag )
#l = plt.legend( loc='best', fontsize=14, frameon=False, handletextpad=0.5 )
#for line, text in zip( l.get_lines(), l.get_texts() ):
#   text.set_color(line.get_color())
fig.savefig( outdir + '/' + figname )
plt.close()



figname = 'CavitySize_vs_t.pdf'
print '\n'+figname
fig = plt.figure( figsize=(5,4), dpi=600 )
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.17, right=right, bottom=bottom, top=top)
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
   item.set_fontsize(fsz_num)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
   item.set_fontsize(fsz_axis)
for axis in ['top','bottom','left','right']:
   ax.spines[axis].set_linewidth(lw_ax)
ax.xaxis.set_tick_params(width=lw_ax)
ax.yaxis.set_tick_params(width=lw_ax)
ax.set_ylabel(r'${\rm log}_{10}(R_{\rm cavity}/{\rm kpc})$');  #ax.set_ylim(ylim)
ax.set_xlabel(r'${\rm Time \/ (Myr)}$');  ax.set_xlim(0,xlim[1])
for stag, scolor in zip( tag_list, sim_col ):
   thesim = sim[stag]
   t = thesim['time'][:] - time_start
   Mgas = thesim['sph']['GasMass'][:,:]
   sigma_gas = np.log10(Mgas/dA)
   indR = np.argmax( sigma_gas, axis=1 )
   #Ngas_inbin = thesim['sph']['Ngas'][:,:]
   #indR = (Ngas_inbin!=0).argmax(axis=1)
   R_cavity = x[indR]
   if stag == 'nof':
     ax.plot( t, R_cavity, '-', color=nof_color, lw=nof_lw, label=stag )
   else:
     ax.plot( t, R_cavity, '-', color=scolor, lw=lw, label=stag )
l = plt.legend( loc='best', fontsize=14, frameon=False, handletextpad=0.5 )
for line, text in zip( l.get_lines(), l.get_texts() ):
   text.set_color(line.get_color())
fig.savefig( outdir + '/' + figname )
plt.close()



figname = 'Rwind_vs_t.pdf'
print '\n'+figname
fig = plt.figure( figsize=(5,4), dpi=600 )
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.17, right=right, bottom=bottom, top=top)
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
   item.set_fontsize(fsz_num)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
   item.set_fontsize(fsz_axis)
for axis in ['top','bottom','left','right']:
   ax.spines[axis].set_linewidth(lw_ax)
ax.xaxis.set_tick_params(width=lw_ax)
ax.yaxis.set_tick_params(width=lw_ax)
ax.set_ylabel(r'${\rm log}_{10}(R_{\rm wind}/{\rm kpc})$');  #ax.set_ylim(ylim)
ax.set_xlabel(r'${\rm Time \/ (Myr)}$');  #ax.set_xlim(0,xlim[1])
for stag, scolor in zip( tag_list, sim_col ):
   if stag == 'nof':
     continue
   thesim = sim[stag]
   t = thesim['time'][:] - time_start
   Ngas_inbin = thesim['sph']['Ngas_BAL'][:,::-1]
   indR = (Ngas_inbin!=0).argmax(axis=1)
   R_wind = x[::-1][indR]
   ax.plot( t, R_wind, '-', color=scolor, lw=lw, label=stag )
   #ax.plot( np.log10(1+t), R_wind, '-', color=scolor, lw=lw, label=stag )
l = plt.legend( loc='best', fontsize=14, frameon=False, handletextpad=0.5 )
for line, text in zip( l.get_lines(), l.get_texts() ):
   text.set_color(line.get_color())
fig.savefig( outdir + '/' + figname )
plt.close()



figname = 'MdotOUT_vs_t.pdf'
print '\n'+figname
fig = plt.figure( figsize=(5,4), dpi=600 )
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.17, right=right, bottom=bottom, top=top)
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
   item.set_fontsize(fsz_num)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
   item.set_fontsize(fsz_axis)
for axis in ['top','bottom','left','right']:
   ax.spines[axis].set_linewidth(lw_ax)
ax.xaxis.set_tick_params(width=lw_ax)
ax.yaxis.set_tick_params(width=lw_ax)
ax.set_ylabel(r'${\rm log}_{10}(\dot{M}_{\rm out}/{\rm M}_{\odot}{\rm yr}^{-1})$');  #ax.set_ylim(ylim)
ax.set_xlabel(r'${\rm Time \/ (Myr)}$');  ax.set_xlim(xlim)
for stag, scolor in zip( tag_list, sim_col ):
   thesim = sim[stag]
   t = thesim['time'][:] - time_start
   rcut = 1.
   indR = np.argmin( np.abs(radial_bins[1:] - rcut) )
   OutflowRate = thesim['sph']['GasOutflowRate'][:,indR] + thesim['sph']['GasInflowRate'][:,indR]
   OutflowRate[OutflowRate<0] = 1e-3
   OutflowRate_BAL = thesim['sph']['GasOutflowRate_BAL'][:,indR]
   #Mout = np.log10( OutflowRate + OutflowRate_BAL )
   #Mout = np.log10( OutflowRate_BAL )
   Mout = np.log10( OutflowRate )
   if stag == 'nof':
     ax.plot( t, Mout, '-', color=nof_color, lw=nof_lw, label=stag )
   else:
     ax.plot( t, Mout, '-', color=scolor, lw=lw, label=stag )
l = plt.legend( loc='best', fontsize=14, frameon=False, handletextpad=0.5 )
for line, text in zip( l.get_lines(), l.get_texts() ):
   text.set_color(line.get_color())
fig.savefig( outdir + '/' + figname )
plt.close()


figname = 'MdotOUT_vs_r.pdf'
print '\n'+figname
fig = plt.figure( figsize=(5,4), dpi=600 )
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.17, right=right, bottom=bottom, top=top)
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
   item.set_fontsize(fsz_num)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
   item.set_fontsize(fsz_axis)
for axis in ['top','bottom','left','right']:
   ax.spines[axis].set_linewidth(lw_ax)
ax.xaxis.set_tick_params(width=lw_ax)
ax.yaxis.set_tick_params(width=lw_ax)
ax.set_ylabel(r'${\rm log}_{10}(\dot{M}_{\rm out}/{\rm M}_{\odot}{\rm yr}^{-1})$');  ax.set_ylim([-1.5,3.5])
ax.set_xlabel(r'${\rm log}_{10}(R \/ / \/ {\rm kpc})$');  ax.set_xlim([-2,0])
k = 1
for stag, scolor in zip( tag_list, sim_col ):
   if stag == 'nof':
     continue
   thesim = sim[stag]
   t = thesim['time'][:] - time_start
   indt = np.where(t>=0)[0]
   #indt = np.where((t>0)&(t<1))[0]

   #OutflowRate = np.mean( thesim['sph']['GasOutflowRate'][indt,:] + thesim['sph']['GasInflowRate'][indt,:], axis=0 )
   #OutflowRate[OutflowRate<0] = 1e-3
   #Mout = np.log10( OutflowRate + OutflowRate_BAL )
   OutflowRate_BAL = np.zeros(nbin)
   GasOutflowVel = deepcopy(thesim['sph']['GasOutflowVel_BAL'][:,:])
 #  GasOutflowVel[ GasOutflowVel>1e4 ] = 1e4
   GasMassOutflow = thesim['sph']['GasMassOutflow_BAL'][:,:]
   for i in range(nbin):
      #indc = np.where( thesim['sph']['GasOutflowVel_BAL'][indt,i] < 5e4 )[0]
      #OutflowRate_BAL[i] = np.mean( thesim['sph']['GasOutflowRate_BAL'][indt[indc],i], axis=0 )
      #OutflowRate_BAL[i] = np.median( thesim['sph']['GasOutflowRate_BAL'][indt,i] )
      dRbin = radial_bins[i+1] - radial_bins[i]
      OutflowRate_BAL[i] = np.mean( GasMassOutflow[indt,i] * (GasOutflowVel[indt,i]*CM_PER_KM/CM_PER_KPC*SEC_PER_YEAR) / dRbin )
   Mout = np.log10( OutflowRate_BAL )
   ax.plot( x, Mout, '-', color=scolor, lw=lw, label=stag )

   mdot = mdot_list[k]
   facc = facc_list[k]   
   Mout_BAL = np.log10( (1.-facc)/facc * mdot )
   ax.plot( [-2,0], [Mout_BAL,Mout_BAL], '--', color=scolor, lw=lw )
   k += 1
l = plt.legend( loc='best', fontsize=14, frameon=False, handletextpad=0.5 )
for line, text in zip( l.get_lines(), l.get_texts() ):
   text.set_color(line.get_color())
fig.savefig( outdir + '/' + figname )
plt.close()


figname = 'VelOUT_vs_r.pdf'
print '\n'+figname
fig = plt.figure( figsize=(5,4), dpi=600 )
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.17, right=right, bottom=bottom, top=top)
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
   item.set_fontsize(fsz_num)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
   item.set_fontsize(fsz_axis)
for axis in ['top','bottom','left','right']:
   ax.spines[axis].set_linewidth(lw_ax)
ax.xaxis.set_tick_params(width=lw_ax)
ax.yaxis.set_tick_params(width=lw_ax)
ax.set_ylabel(r'$V_{\rm wind} ({\rm km\/s}^{-1})$');  ax.set_ylim([2,5.5])
ax.set_xlabel(r'${\rm log}_{10}(R \/ / \/ {\rm kpc})$');  ax.set_xlim([-2,0])
for stag, scolor in zip( tag_list, sim_col ):
   if stag == 'nof':
     continue
   thesim = sim[stag]
   t = thesim['time'][:] - time_start
   indt = np.where(t>=0)[0]
   #indt = np.where((t>0)&(t<1))[0]
   #Vel = np.mean( thesim['sph']['GasOutflowVel'][indt,:] )
   Vel_BAL = np.mean( thesim['sph']['GasOutflowVel_BAL'][indt,:], axis=0 )
   ax.plot( x, np.log10( Vel_BAL ), '-', color=scolor, lw=lw, label=stag )
   #cs_BAL = np.mean( thesim['sph']['SoundSpeed_BAL'][indt,:], axis=0 )
   #ax.plot( x, np.log10( cs_BAL ), '--', color=scolor, lw=lw )
ax.plot( [-2,0], [1e4,1e4], '--', color='black', lw=lw )
l = plt.legend( loc='best', fontsize=14, frameon=False, handletextpad=0.5 )
for line, text in zip( l.get_lines(), l.get_texts() ):
   text.set_color(line.get_color())
fig.savefig( outdir + '/' + figname )
plt.close()



figname = 'Vwind_at_Rwind_vs_t.pdf'
print '\n'+figname
fig = plt.figure( figsize=(5,4), dpi=600 )
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.17, right=right, bottom=bottom, top=top)
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
   item.set_fontsize(fsz_num)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
   item.set_fontsize(fsz_axis)
for axis in ['top','bottom','left','right']:
   ax.spines[axis].set_linewidth(lw_ax)
ax.xaxis.set_tick_params(width=lw_ax)
ax.yaxis.set_tick_params(width=lw_ax)
ax.set_ylabel(r'$V_{\rm wind} ({\rm km\/s}^{-1})$');  #ax.set_ylim(ylim)
ax.set_xlabel(r'${\rm Time \/ (Myr)}$');  #ax.set_xlim(0,xlim[1])
for stag, scolor in zip( tag_list, sim_col ):
   if stag == 'nof':
     continue
   thesim = sim[stag]
   t = thesim['time'][:] - time_start
   Ngas_inbin = thesim['sph']['Ngas_BAL'][:,::-1]
   indR = (Ngas_inbin!=0).argmax(axis=1)
   #R_wind = x[::-1][indR]
   Vel = thesim['sph']['GasOutflowVel_BAL'][:,::-1]
   V_wind = np.zeros(indR.size)
   for i in range(indR.size):
      V_wind[i] = Vel[i,indR[i]]
   V_wind_test = np.diagonal(Vel[:,indR])
   #ax.plot( t, V_wind, '-', color=scolor, lw=lw, label=stag )
   ax.plot( np.log10(1+t), V_wind, '-', color=scolor, lw=lw, label=stag )
l = plt.legend( loc='best', fontsize=14, frameon=False, handletextpad=0.5 )
for line, text in zip( l.get_lines(), l.get_texts() ):
   text.set_color(line.get_color())
fig.savefig( outdir + '/' + figname )
plt.close()




for thistime in [10.]:

 figname = 'SigmaStar_vs_r_%2dMyr.pdf' % int(thistime)
 print '\n'+figname
 fig = plt.figure( figsize=(5,4), dpi=600 )
 ax = fig.add_subplot(111)
 fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
 for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(fsz_num)
 for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
    item.set_fontsize(fsz_axis)
 for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(lw_ax)
 ax.xaxis.set_tick_params(width=lw_ax)
 ax.yaxis.set_tick_params(width=lw_ax)
 ylim = [9,13]
 xlim = [-2,0]
 ax.set_ylabel(r'${\rm log}_{10}(\Sigma_{\rm star}/{\rm M}_{\odot}{\rm kpc}^{-2})$');  ax.set_ylim(ylim);  ax.set_yticks([9,10,11,12,13])
 ax.set_xlabel(r'${\rm log}_{10}(R \/ / \/ {\rm kpc})$');  ax.set_xlim(xlim)
 for stag, scolor in zip( tag_list, sim_col ):
    thesim = sim[stag]
    t = thesim['time'][:] - time_start
    i = np.argmin( np.abs(t - thistime) )
    print 'time =', t[i], 'Myr'
    Mstar = thesim['cyl']['1kpc']['StarMass'][i,:]
    ind = np.where( thesim['cyl']['1kpc']['Nstar'][i,:] > 0 )[0]
    xstar = x[ind]
    sigma_star = np.log10(Mstar[ind]/dA[ind])
    if stag == 'nof':
      ax.plot( xstar, sigma_star, '-', color=nof_color, lw=nof_lw, label=stag )
    else:
      ax.plot( xstar, sigma_star, '-', color=scolor, lw=lw, label=stag )
 l = plt.legend( loc='best', fontsize=14, frameon=False, handletextpad=0.5 )
 for line, text in zip( l.get_lines(), l.get_texts() ):
    text.set_color(line.get_color())
 fig.savefig( outdir + '/' + figname )
 plt.close()

 
 figname = 'SigmaGas_vs_r_%2dMyr.pdf' % int(thistime)
 print '\n'+figname
 fig = plt.figure( figsize=(5,4), dpi=600 )
 ax = fig.add_subplot(111)
 fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
 for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(fsz_num)
 for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
    item.set_fontsize(fsz_axis)
 for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(lw_ax)
 ax.xaxis.set_tick_params(width=lw_ax)
 ax.yaxis.set_tick_params(width=lw_ax)
 ylim = [7,11]
 xlim = [-2,0]
 ax.set_ylabel(r'${\rm log}_{10}(\Sigma_{\rm gas}/{\rm M}_{\odot}{\rm kpc}^{-2})$');  ax.set_ylim(ylim);  ax.set_yticks([7,8,9,10,11])
 ax.set_xlabel(r'${\rm log}_{10}(R \/ / \/ {\rm kpc})$');  ax.set_xlim(xlim)
 for stag, scolor in zip( tag_list, sim_col ):
    thesim = sim[stag]
    t = thesim['time'][:] - time_start
    i = np.argmin( np.abs(t - thistime) )
    #Mgas = thesim['cyl']['1kpc']['GasMass'][i,:]
    #ind = np.where( thesim['cyl']['1kpc']['Ngas'][i,:] > 0 )[0]
    Mgas = thesim['sph']['GasMass'][i,:]
    Ngas_inbin = thesim['sph']['Ngas'][i,:]
    ind = np.where( Ngas_inbin > 0 )[0]
    xgas = x[ind]
    sigma_gas = np.log10(Mgas[ind]/dA[ind])
    if stag == 'nof':
      ax.plot( xgas, sigma_gas, '-', color=nof_color, lw=nof_lw, label=stag )
    else:
      ax.plot( xgas, sigma_gas, '-', color=scolor, lw=lw, label=stag )
 l = plt.legend( loc='best', fontsize=14, frameon=False, handletextpad=0.5 )
 for line, text in zip( l.get_lines(), l.get_texts() ):
    text.set_color(line.get_color())
 fig.savefig( outdir + '/' + figname )
 plt.close()











print '\nDone!!\n'




