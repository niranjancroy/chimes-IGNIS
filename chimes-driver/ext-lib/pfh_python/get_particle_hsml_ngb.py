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
import gadget as gg
import matplotlib.pylab as pylab
import matplotlib.pyplot as plot
import scipy.fftpack as fft
import scipy.misc
import scipy.interpolate as interpolate
import scipy.optimize as optimize
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import ndimage

#
# this is the 'wrapper' routine to the compiled C-library "neighbor_finder" (which should 
#  be compiled with everything else by the shell script "make_all_pylibs.sh"). here you feed 
#  in some list of search center locations, some list of coordinates of the 'target' particles
#  which are eligible to be neighbors, and then optional parameters (weights for the neighbors, 
#  box-sizes for periodic boxes, number of dimensions to project the data set to, minimum/maximum 
#  search lengths for neighbors, target neighbor numbers) -- all the optional parameters have
#  defaults which will be calculated or set for you if you don't set them. the c-routine then 
#  builds a search tree, find neighbors for each particle, and by default returns the kernel-weighted 
#  sum (essentially the weighted number density of neighbor particles) around each search location. 
#  it can of course be modified to calculate any other kernel-weighted quantity. the routine then 
#  returns the search-lengths and sums it found.
#

# these are some random 'utility' routines needed by the master routine below
def snap_ext(snum,four_char=0):
	ext='00'+str(snum);
	if (snum>=10): ext='0'+str(snum)
	if (snum>=100): ext=str(snum)
	if (four_char==1): ext='0'+ext
	if (snum>=1000): ext=str(snum)
	return ext;
def fcor(x):
    return np.array(x,dtype='f',ndmin=1)
def vfloat(x):
    return x.ctypes.data_as(ctypes.POINTER(ctypes.c_float));
def cfloat(x):
    return ctypes.c_float(x);
def checklen(x):
    return len(np.array(x,ndmin=1));
def ok_scan(input,xmax=1.0e10,pos=0):
    if (pos==1):
        return (np.isnan(input)==False) & (abs(input)<=xmax) & (input > 0.);
    else:
        return (np.isnan(input)==False) & (abs(input)<=xmax);

# here is the actual master routine
def get_particle_hsml_ngb( 
        x_search, y_search, z_search, # coordinates of points at which to center the evaluation
        x_target, y_target, z_target, # coordinates of points to evaluate the number density of
        weight_target=0, # weights to assign target points (default is uniform, all equal, but can set to e.g. masses)
        boxsize_x=1.e37, boxsize_y=1.e37, boxsize_z=1.e37, # manually set box side-lengths (for periodic data)
        Numdims=0, # 0=determine from data
        Hmin=0., # default to zero, but can set to some desired value
        Hmax=0., # 0=determine from data
        DesNgb=0 # 0=determine based on dimension (1d=4, 2d=12, 3d=32)
        ):
        
    # project down the data set (should be true already, but in case it isnt)
    if(Numdims==1): 
        y_target*=0; z_target*=0; y_search*=0; z_search*=0;
    if(Numdims==2):
        z_target*=0; z_search*=0;
    # format the data set and remove bad points
    x_search=fcor(x_search); y_search=fcor(y_search); z_search=fcor(z_search); N_search=checklen(x_search); 
    ok=(ok_scan(x_search) & ok_scan(y_search) & ok_scan(z_search)); x_search=x_search[ok]; y_search=y_search[ok]; z_search=z_search[ok];
    x_target=fcor(x_target); y_target=fcor(y_target); z_target=fcor(z_target); N_target=checklen(x_target); 
    ok=(ok_scan(x_target) & ok_scan(y_target) & ok_scan(z_target)); x_target=x_target[ok]; y_target=y_target[ok]; z_target=z_target[ok];
    # calculate the maximum of coordinates still in the box
    if(Hmax==0.):
        dx=np.max(x_target)-np.min(x_target); dy=np.max(y_target)-np.min(y_target); dz=np.max(z_target)-np.min(z_target); 
        Hmax=np.max([dx,dy,dz]); 
    # if numdims not set, guess it automatically from the data
    if(Numdims==0):
        Numdims=3;
        dx=np.max(x_target)-np.min(x_target); 
        if(dx<=0): Numdims-=1;
        dy=np.max(y_target)-np.min(y_target); 
        if(dy<=0): Numdims-=1;
        dz=np.max(z_target)-np.min(z_target); 
        if(dz<=0): Numdims-=1;
    # if desngb not set, resort to defaults
    if(DesNgb==0):
        if(Numdims==1): DesNgb=4;
        if(Numdims==2): DesNgb=8; #12;
        if(Numdims==2): DesNgb=16; #12;
        if(Numdims==3): DesNgb=16;#32;
        if(Numdims==3): DesNgb=64;#32;
    # if weights are not set, default to uniform weights
    wt=weight_target
    if(checklen(weight_target)<=1):
        wt = 0.*x_target + 1.;
    else:
        wt = fcor(wt); 
        wt = wt[ok];

    ## load the c-routine we need, and prepare data casts
    exec_call=util.return_python_routines_cdir()+'/neighbor_finder/kde.so'
    h_routine=ctypes.cdll[exec_call];
    h_out_cast=ctypes.c_float*N_search; H_OUT=h_out_cast(); NGB_OUT=h_out_cast(); DesNgb=np.int(DesNgb);
    ## main call to the hsml-finding routine
    h_routine.kde( ctypes.c_int(N_search), vfloat(x_search), vfloat(y_search), vfloat(z_search), \
        ctypes.c_int(N_target), vfloat(x_target), vfloat(y_target), vfloat(z_target), vfloat(wt), \
        cfloat(boxsize_x), cfloat(boxsize_y), cfloat(boxsize_z), \
        ctypes.c_int(Numdims), ctypes.c_int(DesNgb), cfloat(Hmin), cfloat(Hmax), \
        ctypes.byref(H_OUT), ctypes.byref(NGB_OUT) )

    ## now put the output arrays into a useful format 
    h_smoothing = np.double(np.ctypeslib.as_array(np.copy(H_OUT)));
    numden_neighbors = np.double(np.ctypeslib.as_array(np.copy(NGB_OUT)));
    ## normalize quantities to the box-scale lengths
    Hmax = np.double(Hmax); N_target = np.double(N_target);
    h_smoothing /= boxsize_x; # Normalize to x length, this is taken from header in grain_density_from_snapshot.py
    #  JS: Using weights=masses means we don't need to normalize by anything to get the mass density. If you want number density instead, uncomment this.
    if(checklen(weight_target)<=1):
        numden_neighbors /= (N_target/(boxsize_x*boxsize_y*boxsize_z));
    
    return h_smoothing, numden_neighbors, Numdims, DesNgb

