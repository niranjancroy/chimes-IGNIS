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
import gadget
import matplotlib.pylab as pylab
import matplotlib.pyplot as plot
import scipy.fftpack as fft
import scipy.misc
import scipy.interpolate as interpolate
import scipy.optimize as optimize
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

def test(snum=60):
    rmin = 1.e-5
    pylab.close('all')
    pylab.yscale('log'); pylab.xscale('linear'); pylab.axis([np.log10(rmin),-0.4,1.e-2,100.])

    for ptype,color,linesty in zip([0,3],['blue','red'],['--','-']):
        P=gadget.readsnap('/Users/phopkins/Downloads/grains',snum,ptype)
        N0 = P['p'][:,0].size
        ok = np.where(np.random.rand(N0) < 1.e4/N0)
        r,n = get_particle_corrfn(P['p'][:,0][ok],P['p'][:,1][ok],P['p'][:,2][ok],
            rmin,100.)
        print('n == ',n)
        pylab.plot(r,n,color=color,linewidth=2.,linestyle=linesty)
        # 41 (2D, alpha=0.01) - dust pwr-law to ~1e-4 and below, slope very close to xi~r^(-1/2), gas flattens at 1e-2
        # 60 (3D, alpha=0.01) - similar pwr-law (maybe a bit more shallow; difference smaller with gas)
        # 120 - same run as 60; difference now very clear. 
        # 57 (3d, alpha=0.1) - similar to 0.01, very clear
        # 55 (3d, alpha=1) - much lower power, but similar power-law 

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



# here is the routine which synthesizes from the master routine
def get_particle_corrfn( 
        x_search, y_search, z_search, # coordinates of points at which to center the evaluation
        rmin, N_rgrid):

    N_search_min = 0
    N_search_max = 1000000
    r_grid_out, n_grid_out = get_particle_corrfn_workhorse(x_search,y_search,z_search,
        N_search_min,N_search_max,rmin,N_rgrid)
    
    N_search = x_search.size;    
    N_total = N_search*(N_search-1)/2.;
    n_grid_out /= N_total;
    
    ## now we need to turn it into an actual correlation function (currently just 
    ##   number per log interval in distance)
    numdims = 3
    if(np.max(z_search)-np.min(z_search)==0): numdims = 2
    print('Ndims = ',numdims)
    dVol=0.*r_grid_out
    rg=10.**r_grid_out
    if(numdims==3): Vol=4.*np.pi/3.*rg*rg*rg
    if(numdims==2): Vol=np.pi*rg*rg
    dVol[0]=Vol[0]; dVol[1::]=np.diff(Vol);
    n_normalized = n_grid_out / dVol - 1
    log_r_grid = 1.*r_grid_out

    return log_r_grid, n_normalized


def dump_corrfn(sdir='/Users/phopkins/Downloads/grains/',snum=57,
        Nmin=0, Nmax=100000, rmin=1.e-5, N_rgrid=100, ptype=3):
    P=gadget.readsnap(sdir,snum,ptype)
    N0 = P['p'][:,0].size
    ok = np.where(np.random.rand(N0) < 1.e5/N0)
    s=np.argsort(np.random.rand(N0)); ok=s[ok];
    x=P['p'][:,0][ok]; y=P['p'][:,1][ok]; z=P['p'][:,2][ok];
    if((2==2)):#&(ptype==0)):
        n_subsample = 10;
        h_eff = 3.0 * (P['m'][ok]/P['rho'][ok])**(1./3.)
        #h_eff = 1.e-3
        xx=np.zeros(0); yy=np.zeros(0); zz=np.zeros(0); n0=x.size;
        for i in range(n_subsample):
            xx=np.append(xx,x+h_eff*np.random.normal(0.,1.,n0))
            yy=np.append(yy,y+h_eff*np.random.normal(0.,1.,n0))
            zz=np.append(zz,z+h_eff*np.random.normal(0.,1.,n0))
            #xx=np.append(xx,x+h_eff*(np.random.rand(n0)-0.5))
            #yy=np.append(yy,y+h_eff*(np.random.rand(n0)-0.5))
            #zz=np.append(zz,z+h_eff*(np.random.rand(n0)-0.5))
        x=xx; y=yy; z=zz;
    r,n = get_particle_corrfn_workhorse(x,y,z,
        Nmin, Nmax, rmin, N_rgrid, snapnum=snum, saveoutput=True, sdir=sdir)


# here is the actual master routine
def get_particle_corrfn_workhorse( 
        x_search, y_search, z_search, # coordinates of points at which to center the evaluation
        N_search_min, N_search_max, rmin, N_rgrid, 
        snapnum=0, saveoutput=True, sdir='/Users/phopkins/Downloads/grains/'):
        
    # format the data set and remove bad points
    x_search=fcor(x_search); y_search=fcor(y_search); z_search=fcor(z_search); N_search=checklen(x_search); 
    ok=(ok_scan(x_search) & ok_scan(y_search) & ok_scan(z_search)); x_search=x_search[ok]; y_search=y_search[ok]; z_search=z_search[ok];

    ## load the c-routine we need, and prepare data casts
    exec_call=util.return_python_routines_cdir()+'/CorrFn/corrfn.so'
    h_routine=ctypes.cdll[exec_call];
    N_rgrid=np.int(N_rgrid); N_search_min=np.int(N_search_min); N_search_max=np.int(N_search_max);
    h_out_cast=ctypes.c_float*N_rgrid; R_OUT=h_out_cast(); N_OUT=h_out_cast(); 
    ## main call to the hsml-finding routine
    h_routine.kde( ctypes.c_int(N_search), vfloat(x_search), vfloat(y_search), vfloat(z_search), \
        ctypes.c_int(N_search_min), ctypes.c_int(N_search_max), \
        cfloat(rmin), ctypes.c_int(N_rgrid), \
        ctypes.byref(R_OUT), ctypes.byref(N_OUT) )

    ## now put the output arrays into a useful format 
    r_grid_out = np.double(np.ctypeslib.as_array(np.copy(R_OUT)));
    n_grid_out = np.double(np.ctypeslib.as_array(np.copy(N_OUT)));

    if(saveoutput):
        snap_string=snap_ext(snapnum,four_char=1)
        nmin_string=snap_ext(N_search_min/1000,four_char=1)
        nmax_string=snap_ext(N_search_max/1000,four_char=1)

        filename=sdir+'/tmp_corrfn_'+snap_string+'_Nmin_'+nmin_string+'_Nmax_'+nmax_string+'.dat'
        outfi = h5py.File(filename,'w')
        dset_xvec = outfi.create_dataset('rgrid',data=r_grid_out)
        dset_yvec = outfi.create_dataset('ngrid',data=n_grid_out)
        outfi.close()            


    return r_grid_out, n_grid_out

