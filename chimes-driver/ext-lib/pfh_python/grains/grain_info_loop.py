import numpy as np
import ctypes
import pfh_utils as util
import math
import os.path
import struct
import array
import h5py as h5py
import os.path
import math
import matplotlib
#matplotlib.use('Agg') ## this calls matplotlib without a X-windows GUI
import gadget as gadget
import matplotlib.pylab as pylab
import matplotlib.pyplot as plot
import scipy.fftpack as fft
import scipy.misc
import scipy.interpolate as interpolate
import scipy.optimize as optimize
import get_particle_hsml_ngb
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import ndimage


def grain_info_loop():
    s0='/scratch/03616/tg829152/'
    compile_grain_info_loop(sdir=s0+'rhospec_mach10',snum_min=15,snum_max=16)


def compile_grain_info_loop(sdir='./',snum_min=1,snum_max=1):
    hmin_vec = np.array([0.,1.e-3,1.e-2])
    xvec = np.arange(-8.,6.,0.05); y_total = np.zeros((xvec.size-1,xvec.size-1)); n_total = 0.0;
    # loop over snapshots, and for each, loop over hmin values, to dump the results we want
    for hmin in hmin_vec:
        for snum in np.arange(snum_min,snum_max+1):
            y, xedge, yedge = compile_info_snapshots_sub(xvec,sdir=sdir,snum=snum,Hmin=hmin);
            y_total += y;
            n_total += 1.0;
    y_total /= n_total;
    print(y_total)
    s0=sdir.split("/"); snapdir_specific=s0[len(s0)-1]; n_s0=1;
    if(len(snapdir_specific)<=1): snapdir_specific=s0[len(s0)-2]; 
    n_s0=2; ss=snap_ext(snum,four_char=1)
    hmin_string = "%0.0e" % hmin
    s_out = '/Users/phopkins/Documents/work/plots/grain_turb_box/rhodust_rhogas_compiled/'
    filename=s_out+snapdir_specific+'_rho_compiled_'+hmin_string+'.dat'
    print(filename)
    outfi = h5py.File(filename,'w')
    dset_xvec = outfi.create_dataset('Histogram_Bins',data=xvec)
    dset_yvec = outfi.create_dataset('Histogram_Vals',data=y_total)
    outfi.close()            
    return xvec, y_total


def compile_info_snapshots_sub(xvec,sdir='grains',snum=2,Hmin=0.0):
    s0=sdir.split("/"); snapdir_specific=s0[len(s0)-1]; n_s0=1;
    if(len(snapdir_specific)<=1): snapdir_specific=s0[len(s0)-2]; n_s0=2;
    ss=snap_ext(snum,four_char=1)
    hmin_string = "%0.0e" % Hmin
    filename=sdir+'/'+snapdir_specific+'_ngbfile_hmin_'+hmin_string+'_snap_'+ss+'.dat'
    print(filename)
    outfi = h5py.File(filename,'r')
    nngb_d    = np.array(outfi['Mass_Density_List_Of_Dust_Neighbors'])
    nngb_g    = np.array(outfi['Mass_Density_List_Of_Gas_Neighbors'])
    outfi.close()
    ok=np.where((np.isnan(nngb_g)==False)&(np.isnan(nngb_d)==False)&(nngb_g>0)&(nngb_d>0))
    #nngb_d[ok] *= (8./np.pi) / (4./3.) # correction for bug in old code version
    xin = np.log10(nngb_g[ok])
    yin = np.log10(nngb_d[ok])
    info_hist_2d, xedges, yedges = np.histogram2d(xin, yin, bins=(xvec,xvec), normed=True)

    return info_hist_2d, xedges, yedges


# quick subroutine to convert snapshot number to a file extension string
def snap_ext(snum,four_char=0):
	ext='00'+str(snum);
	if (snum>=10): ext='0'+str(snum)
	if (snum>=100): ext=str(snum)
	if (four_char==1): ext='0'+ext
	if (snum>=1000): ext=str(snum)
	return ext;

