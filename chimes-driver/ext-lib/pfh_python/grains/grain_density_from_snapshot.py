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
    s0='./'; snum=303; hmin=0.
    #plot_grain_density_from_snapshot(snum=303)
    xvec = np.arange(-8.,6.,0.05); y_total = 0.0*xvec; 
    y, xedge, yedge = compile_info_snapshots_sub(xvec,sdir=s0,snum=snum,Hmin=hmin)


def compile_grain_info_loop(sdir='./',snum_min=1,snum_max=1):
    hmin_vec = np.array([0.])
    xvec = np.arange(-8.,6.,0.05); y_total = 0.0*xvec; 
    # loop over snapshots, and for each, loop over hmin values, to dump the results we want
    for hmin in hmin_vec:
        for snum in np.arange(snum_min,snum_max+1):
            y, xedge, yedge = compile_info_snapshots_sub(xvec,sdir=sdir,snum=snum,Hmin=hmin)
            y_total += y
            
        s0=sdir.split("/"); snapdir_specific=s0[len(s0)-1]; n_s0=1;
        if(len(snapdir_specific)<=1): snapdir_specific=s0[len(s0)-2]; n_s0=2;
        ss=snap_ext(snum,four_char=1)
        hmin_string = "%0.0e" % Hmin
        filename='./'+snapdir_specific+'_rho_compiled_'+hmin_string+'.dat'
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
    filename=snapdir_specific+'_ngbfile_hmin_'+hmin_string+'_snap_'+ss+'.dat'
    #filename='grains_ngbfile_hmin_0e+00_a1e-6_1000_snap_0042.dat'
    print(filename)
    outfi = h5py.File(filename,'r')
    nngb_d    = np.array(outfi['Mass_Density_List_Of_Dust_Neighbors'])
    nngb_g    = np.array(outfi['Mass_Density_List_Of_Gas_Neighbors'])
    h0 = np.array(outfi['Smoothing_Length_List_Of_Dust_Neighbors'])
    outfi.close()
    ok=np.where((np.isnan(nngb_g)==False)&(np.isnan(nngb_d)==False)&(nngb_g>0)&(nngb_d>0))
    bad = np.where((np.isnan(nngb_g)==True)|(np.isnan(nngb_d)==True)|(nngb_d<=0)|(nngb_d<=0))
    print(nngb_d[bad])
    print(nngb_g[bad])
    print(h0[bad])
    nngb_d[bad] = 1.0e-12
    #nngb_d[ok] *= (8./np.pi) / (4./3.) # correction for bug in old code version
    print(np.mean(nngb_d), np.mean(nngb_g), nngb_d.size, nngb_g.size, nngb_d[ok].size, np.min(nngb_d), np.min(nngb_d[ok]))
    xin = np.log10(nngb_g[ok])
    yin = np.log10(nngb_d[ok])
    info_hist_2d, xedges, yedges = np.histogram2d(xin, yin, bins=(xvec,xvec), normed=True)

    return info_hist_2d, xedges, yedges


def test_for_kernel_bias(snum=303,hmin=0.0):
    xvec=np.arange(-8.,6.,0.05); 
    y,xe,ye=compile_info_snapshots_sub(xvec,sdir='grains',snum=snum,Hmin=hmin); 
    n0=xe.size; dx=xe[1]-xe[0]; x0=xe[0:n0-1]+dx; xv,yv=np.meshgrid(x0,x0); 
    dn = y*dx*dx
    print(np.sum(dn), np.sum(dn*10.**(-yv)),np.sum(dn*10.**(-yv+xv)))
    print(np.sum(dn), np.sum(dn*10.**(-xv)),np.sum(dn*10.**(-xv+yv)))
    # each of these should sum to 1, if there are no kernel biases, etc
    # [also, for 'y' from compile_info_snapshots_sub, x-axis is dust, y-axis is gas!]



def grain_density_loop():
    s0='/panfs/ds06/sxs/phopkins/GIZMO_tests/grains/densities'
    compile_grain_density_loop(sdir=s0+'2d_64_rho0pt01_mach5',snum_min=1,snum_max=20)


def compile_grain_density_loop(sdir='/Users/phopkins/Downloads/grains',snum_min=1,snum_max=1):
    hmin_vec = np.array([0.])
    xvec = np.arange(-8.,10.,0.1); y_total = 0.0*xvec; 
    # loop over snapshots, and for each, loop over hmin values, to dump the results we want
    for hmin in hmin_vec:
        for snum in np.arange(snum_min,snum_max+1):
            y = compile_density_snapshots_sub(xvec,sdir=sdir,snum=snum,Hmin=hmin)
            y_total += y
            
        s0=sdir.split("/"); snapdir_specific=s0[len(s0)-1]; n_s0=1;
        if(len(snapdir_specific)<=1): snapdir_specific=s0[len(s0)-2]; n_s0=2;
        ss=snap_ext(snum,four_char=1)
        hmin_string = "%0.0e" % Hmin
        filename='./'+snapdir_specific+'_rho_compiled_'+hmin_string+'.dat'
        outfi = h5py.File(filename,'w')
        dset_xvec = outfi.create_dataset('Histogram_Bins',data=xvec)
        dset_xvec = outfi.create_dataset('Histogram_Vals',data=y_total)
        outfi.close()            
    return xvec, y_total

def compile_density_snapshots_sub(xvec,sdir='grains',snum=2,Hmin=0.0):
    s0=sdir.split("/"); snapdir_specific=s0[len(s0)-1]; n_s0=1;
    if(len(snapdir_specific)<=1): snapdir_specific=s0[len(s0)-2]; n_s0=2;
    ss=snap_ext(snum,four_char=1)
    hmin_string = "%0.0e" % Hmin
    filename='./'+snapdir_specific+'_ngbfile_hmin_'+hmin_string+'_snap_'+ss+'.dat'
    print(filename)
    outfi = h5py.File(filename,'r')
    SnapDir   = outfi["Snapshot_Directory"].value
    SnapNum   = outfi['Snapshot_Number'].value
    SimTime   = outfi['Snapshot_Time'].value
    NumDims   = outfi['Number_Of_Simulation_Dimensions'].value
    BoxSize   = outfi['Periodic_BoxSize'].value
    NdustTot  = outfi['Total_Number_Of_Dust_Particles'].value
    NgasTot   = outfi['Total_Number_Of_Gas_Particles'].value
    Hmin_calc = outfi['Minimum_Smoothing_Length_Allowed'].value
    desnumngb = outfi['Desired_Number_Of_Neighbors'].value
    h_d       = np.array(outfi['Smoothing_Length_List_Of_Dust_Neighbors'])
    h_g       = np.array(outfi['Smoothing_Length_List_Of_Gas_Neighbors'])
    nngb_d    = np.array(outfi['Mass_Density_List_Of_Dust_Neighbors'])
    nngb_g    = np.array(outfi['Mass_Density_List_Of_Gas_Neighbors'])
    outfi.close()
    ok=np.where((np.isnan(nngb_g)==False)&(np.isnan(nngb_d)==False)&(nngb_g>0)&(nngb_d>0)&(h_d>0)&(h_g>0))
    dust_to_gas = np.log10(nngb_d[ok] / (1.e-30 + nngb_g[ok]))
    y = np.histogram(dust_to_gas, bins=xvec, density=True)
    return y



# this routine loops over snapshots to save the density fields (of dust and gas) to binary files, 
#  for a set of minimum smoothing lengths
def grain_density_from_snapshot_loop(sdir='/Users/phopkins/Downloads/grains',snum_min=1,snum_max=1):
    hmin_vec = np.array([0.])#,1.e-3,1.e-2])
    # loop over snapshots, and for each, loop over hmin values, to dump the results we want
    Pd=gadget.readsnap(sdir,snum_min,3);
    a_vec = np.unique(np.round(Pd['GrainSize']*1.0e6))
    for snum in np.arange(snum_min,snum_max+1):
        print(snum)
        for hmin in hmin_vec:
            for a_grain in a_vec:
                hd,nd,hg,ng=grain_density_from_snapshot(sdir=sdir,snum=snum,Hmin=hmin,DEBUG=False,RANDOM=False,a_grain=a_grain)

# this routine makes a heat-map plot of the results from the density field calculation,
#   from the binary files saved by routines like the one above
def plot_grain_density_from_snapshot(sdir='grains',snum=2,Hmin=0.0,figname='grain_tst.pdf'):
    pylab.close('all');
    plot.figure(1,figsize=(4.,4.)); f00=12
    matplotlib.rcParams.update({'font.size':f00})    

    s0=sdir.split("/"); snapdir_specific=s0[len(s0)-1]; n_s0=1;
    if(len(snapdir_specific)<=1): snapdir_specific=s0[len(s0)-2]; n_s0=2;
    ss=snap_ext(snum,four_char=1)
    hmin_string = "%0.0e" % Hmin
    filename='./'+snapdir_specific+'_ngbfile_hmin_'+hmin_string+'_snap_'+ss+'.dat'
    #filename='grains_ngbfile_hmin_0e+00_a1e-6_1000_snap_0042.dat'
    outfi = h5py.File(filename,'r')
    SnapDir   = outfi["Snapshot_Directory"].value
    SnapNum   = outfi['Snapshot_Number'].value
    SimTime   = outfi['Snapshot_Time'].value
    NumDims   = outfi['Number_Of_Simulation_Dimensions'].value
    BoxSize   = outfi['Periodic_BoxSize'].value
    NdustTot  = outfi['Total_Number_Of_Dust_Particles'].value
    NgasTot   = outfi['Total_Number_Of_Gas_Particles'].value
    Hmin_calc = outfi['Minimum_Smoothing_Length_Allowed'].value
    desnumngb = outfi['Desired_Number_Of_Neighbors'].value
    h_d       = np.array(outfi['Smoothing_Length_List_Of_Dust_Neighbors'])
    h_g       = np.array(outfi['Smoothing_Length_List_Of_Gas_Neighbors'])
    nngb_d    = np.array(outfi['Mass_Density_List_Of_Dust_Neighbors'])
    nngb_g    = np.array(outfi['Mass_Density_List_Of_Gas_Neighbors'])
    outfi.close()
    print(SnapDir,SnapNum,SimTime,NumDims,BoxSize,NdustTot,NgasTot,Hmin_calc,desnumngb)
    dg = nngb_d/nngb_g
    print('min/med/mean/std/max-n_d=',np.min(nngb_d),np.median(nngb_d),np.mean(nngb_d),np.std(nngb_d),np.max(nngb_d),' std-log10-n=',np.std(np.log10(nngb_d)))
    print('min/med/mean/std/max-n_g=',np.min(nngb_g),np.median(nngb_g),np.mean(nngb_g),np.std(nngb_g),np.max(nngb_g),' std-log10-n=',np.std(np.log10(nngb_g)))
    print('min/med/mean/std/max-d/g=',np.min(dg),np.median(dg),np.mean(dg),np.std(dg),np.max(dg),' std-log10-n=',np.std(np.log10(dg)))
    print('min/max h_d =',np.min(h_d),np.max(h_d))
    print('min/max h_g =',np.min(h_g),np.max(h_g))
    ok=np.where((np.isnan(nngb_g)==False)&(np.isnan(nngb_d)==False)&(nngb_g>0)&(nngb_d>0)&(h_d>0)&(h_g>0))

    # scatter plot: simple, but very hard to read if there are a lot of points
    #pylab.plot(np.log10(nngb_g[ok]),np.log10(nngb_d[ok]),marker=',',linestyle='',rasterized=True);
    # better to plot a 'heat map' of points. do this with hexbin. since we care about outliers, 
    #   we should use bin='log' to track the points in the outskirts
    pylab.hexbin(np.log10(nngb_g[ok]),np.log10(nngb_d[ok]),gridsize=60,cmap='Paired',bins='log')
    pylab.xlabel(r'$\log_{10}[\ n_{\rm gas} / \langle n_{\rm gas} \rangle\ ]$',fontsize=f00)
    pylab.ylabel(r'$\log_{10}[\ n_{\rm dust} / \langle n_{\rm dust} \rangle\ ]$',fontsize=f00)
    pylab.subplots_adjust(left=0.05,bottom=0.05,right=1-0.05,top=1-0.05,hspace=0.1,wspace=0.1);
    #pylab.savefig(figname,transparent=True,bbox_inches='tight')    
    

# this routine makes and saves plots of the turbulent box at one instant, calling the 
#   'workhorse' routine below 
def visualize_box_plot(snum=0,zslice=0.0,
        sdir='/Users/phopkins/Downloads/grains/',
        figname='grain_vis.pdf',gas_plot_key='density'):
    pylab.close('all');
    plot.figure(1,figsize=(12.,12.)); f00=12
    matplotlib.rcParams.update({'font.size':f00})
    visualize_turbulent_grain_box(sdir,snum,key=gas_plot_key,zslice=zslice)
    pylab.subplots_adjust(left=0.01,bottom=0.01,right=1-0.01,top=1-0.01,hspace=0.1,wspace=0.1);
    #pylab.savefig(figname,transparent=True,bbox_inches='tight')
    
    
# this routine actually constructs the visualization of the turbulent gas+grain box. 
#   gas is visualized as a color-mapped field. default parameter to map is log(density)
#   the dust is overlaid in points. this is done in slices at z-position 'zslice' (for 2D data
#   this parameter is ignored)
def visualize_turbulent_grain_box(sdir,snum,key='density',zslice=0.0):
    P=gadget.readsnap(sdir,snum,0); Pg=P;
    unit_L = 1 #3.086e19
    unit_M = 1 #9.816e35
    unit_v = 1 #35207.0
    unit_E = unit_v*unit_v*unit_M
    unit_vol = unit_L*unit_L*unit_L
    #if(np.max(P['p'][:,2]-np.min(P['p'][:,2])) <= 0.): unit_vol = unit_L*unit_L
    unit_rho = unit_M/unit_vol
    
    x=P['p'][:,0]; y=P['p'][:,1]; z=P['p'][:,2]; d=P['rho']; h=P['Hsml']
    Bmag = np.sqrt(Pg['B'][:,0]**2+Pg['B'][:,1]**2+Pg['B'][:,2]**2)
    L_part = ((3.*Pg['m']/(4.*np.pi*Pg['rho']))**(1./3.))
    divB = np.log10(1.e-10+np.abs(L_part*Pg['h']*Pg['divB']/Bmag))
    d_vx=np.std(Pg['v'][:,0]); d_vy=np.std(Pg['v'][:,1]); d_vz=np.std(Pg['v'][:,2]); 
    print('B=',np.min(Bmag),np.median(Bmag),np.mean(Bmag),np.max(Bmag))
    print('divB=',np.min(divB),np.median(divB),np.mean(divB),np.max(divB))
    print('rho=',np.min(np.log10(d)),np.median(np.log10(d)),np.mean(np.log10(d)),np.max(np.log10(d)),np.std(np.log10(d)))
    print('rho=',np.min(d),np.median(d),np.mean(d),np.max(d),np.std(d))
    print('B_egy=',np.sqrt(np.mean(Bmag*Bmag*P['m']/P['rho'])/np.mean(P['m']/P['rho'])), np.sum(Bmag*Bmag*P['m']/P['rho'])/(8.*np.pi))
    v2 = d_vx**2+d_vy**2+d_vz**2
    print('v_egy=',d_vx,d_vy,d_vz,np.sqrt(v2), 0.5*np.sum(P['m']*(Pg['v'][:,0]**2+Pg['v'][:,1]**2+Pg['v'][:,2]**2)))
    e_turb = 0.5*np.sum(P['m']*(Pg['v'][:,0]**2+Pg['v'][:,1]**2+Pg['v'][:,2]**2))
    e_therm= np.sum(P['m']*P['u']) * (1.001-1.)
    e_mag  = (1./(8.*np.pi))*np.sum(Bmag*Bmag * P['m']/P['rho']) #* 4.*np.pi
    print('energy: turb / thermal / mag = ',e_turb,e_therm,e_mag)
    print('energy: turb / thermal / mag = ',unit_E*e_turb,unit_E*e_therm,unit_vol*e_mag)
    xmax=np.max([np.max(x),np.max(y),np.max(z)])
    x/=xmax; y/=xmax; z/=xmax;
    x[(x>1)]-=1; x[(x<0)]+=1;
    y[(y>1)]-=1; y[(y<0)]+=1;
    z[(z>1)]-=1; z[(z<0)]+=1;
    if(np.max(z)<=zslice): zslice=np.max(z);
    z[(z-zslice)>0.5]-=1; z[(z-zslice)<-0.5]+=1;
    z-=zslice;
    yg,xg=np.mgrid[0:1:512j, 0:1:512j]
    #yg,xg=np.mgrid[0:1:4096j, 0:1:4096j]
    ok=np.where(np.abs(z) < 0.5*h)
    h_median = np.median(h[ok])
    u = np.log10(d); vmin=0.8*np.min(u[ok]); vmax=0.8*np.max(u[ok])
    print('minmax == ',vmin,vmax)

    Pd=gadget.readsnap(sdir,snum,3); 
    xd=Pd['p'][:,0]; yd=Pd['p'][:,1]; zd=Pd['p'][:,2];
    xd/=xmax; yd/=xmax; zd/=xmax;
    xd[(xd>1)]-=1; xd[(xd<0)]+=1;
    yd[(yd>1)]-=1; yd[(yd<0)]+=1;
    zd[(zd>1)]-=1; zd[(zd<0)]+=1;
    zd[(zd-zslice)>0.5]-=1; zd[(zd-zslice)<-0.5]+=1;
    zd-=zslice;
    okd = np.where(np.abs(zd) < 0.5*h_median)
    
    print('size == ',xd[okd].size)

    if(key=='density'): 
        u = np.log10(d); #vmin=-1.; vmax=1.;
        #vmin=np.min(u); vmax=np.max(u); vmin=-1.; vmax=1.8;
    if(key=='divb'): 
        u = divB; vmin=-4.; vmax=0.;
    if(key=='bfield'): 
        u = Bmag; vmin=np.min(Bmag[ok]); vmax=np.max(Bmag[ok]);
    cmap='jet'
    #cmap='gnuplot2'
    dg=interpolate.griddata((x[ok],y[ok]),u[ok],(xg,yg),method='linear',fill_value=np.median(u[ok]));
    im=pylab.imshow(dg,interpolation='bicubic',vmin=vmin,vmax=vmax,cmap=cmap,extent=(0,1,0,1));
    pylab.plot(xd[okd],yd[okd],marker=',',linestyle='',color='black');#,rasterized=True);
    pylab.show()

    frame1=plot.gca() #get coordinate axes, to use for tick manipulation below
    frame1.axes.xaxis.set_ticklabels([]) #no tick names
    frame1.axes.yaxis.set_ticklabels([]) #no tick names
    frame1.axes.get_xaxis().set_ticks([]) #no ticks
    frame1.axes.get_yaxis().set_ticks([]) #no ticks
    if(2==0):
        ax_ins = inset_axes(frame1,width="100%",height="100%",loc=1,
            bbox_to_anchor=(0.05,0.01,0.90,0.05),bbox_transform=frame1.transAxes,borderpad=0)
        if(key=='divb'):
            cb=pylab.colorbar(im,cax=ax_ins,orientation='horizontal',ticks=[-4,-3,-2,-1,0])
            cb.ax.set_xticklabels(['-4','-3','-2','-1','0'],color='white',fontsize=9)
            cb.ax.tick_params(length=3.5,color='black',pad=-15.)
        else:
            cb=pylab.colorbar(im,cax=ax_ins,orientation='horizontal',ticks=[1.,1.5,2.,2.5])
            cb.ax.set_xticklabels(['$1$','$1.5$','$2$','$2.5$'],color='white',fontsize=10)
            cb.ax.tick_params(length=3.5,color='black',pad=-15.)
    return avvz



# quick subroutine to convert snapshot number to a file extension string
def snap_ext(snum,four_char=0):
	ext='00'+str(snum);
	if (snum>=10): ext='0'+str(snum)
	if (snum>=100): ext=str(snum)
	if (four_char==1): ext='0'+ext
	if (snum>=1000): ext=str(snum)
	return ext;


# this is the 'master' routine to take a snapshot of a turbulent grain box, pull out the relevant 
#   numbers (dust and gas positions), feed it to subroutine get_particle_hsml_ngb which wraps to 
#   the c-libraries that do the neighbor calculations, then saves the outputs to a binary file so 
#   we don't have to re-do the calculation again
def grain_density_from_snapshot(sdir='/Users/phopkins/Downloads/grains',snum=1,Hmin=0.,boxsize_x=1,boxsize_y=1,boxsize_z=1,DEBUG=False,RANDOM=False,a_grain=1.):
    # read in the snapshot data
    Pd=gadget.readsnap(sdir,snum,3); Pg=gadget.readsnap(sdir,snum,0); Ph=gadget.readsnap(sdir,snum,0,header_only=1);
    BoxSize=Ph['boxsize']; time=Ph['time'];
    boxsize_x*=BoxSize; boxsize_y*=BoxSize; boxsize_z*=BoxSize;
    print(boxsize_y)
    if(RANDOM==True):
        ndust=np.size(Pd['p'][:,0]); 
        xmin=np.min(Pd['p'][:,0]); xmax=np.max(Pd['p'][:,0]);
        Pd['p'][:,0]=xmin+(xmax-xmin)*np.random.rand(ndust) # randomize positions
        Pd['p'][:,1]=xmin+(xmax-xmin)*np.random.rand(ndust) # randomize positions
        Pd['p'][:,2]=xmin+(xmax-xmin)*np.random.rand(ndust) # randomize positions
    # center searching on dust locations
    # CHANGE to be grid
    xs=Pd['p'][:,0]; ys=Pd['p'][:,1]; zs=Pd['p'][:,2];
    Ndust = Pd['p'][:,0].size
    Ngas = Pg['p'][:,0].size

    if(DEBUG==True):
        b0=0.2
        ok_d=np.where((Pd['p'][:,0]<b0)&(Pd['p'][:,1]<b0)&(Pd['p'][:,2]<b0))
        ok_g=np.where((Pg['p'][:,0]<b0)&(Pg['p'][:,1]<b0)&(Pg['p'][:,2]<b0))
    else:
        ok_d=np.where(np.abs(Pd['m']) > -1)
        # ok_d=np.where(np.round(Pd['GrainSize']*1.0e6)==a_grain);
        ok_g=np.where(np.abs(Pg['m']) > -1)
    
    DesNgb=0 # Default
    h_d,nngb_d,ndims,Nngb = get_particle_hsml_ngb.get_particle_hsml_ngb(\
                                    xs[ok_d],ys[ok_d],zs[ok_d],Pd['p'][:,0][ok_d],Pd['p'][:,1][ok_d],Pd['p'][:,2][ok_d],\
                                                                        Hmin=Hmin,boxsize_x=boxsize_x,boxsize_y=boxsize_y,boxsize_z=boxsize_z,DesNgb=DesNgb,weight_target=Pd['m'][ok_d])
    h_g,nngb_g,ndims,Nngb = get_particle_hsml_ngb.get_particle_hsml_ngb(\
                                    xs[ok_d],ys[ok_d],zs[ok_d],Pg['p'][:,0][ok_g],Pg['p'][:,1][ok_g],Pg['p'][:,2][ok_g],\
                                                                        Hmin=Hmin,boxsize_x=boxsize_x,boxsize_y=boxsize_y,boxsize_z=boxsize_z,DesNgb=DesNgb,weight_target=Pg['m'][ok_g])
    dg = nngb_d/nngb_g
    print('min/med/mean/std/max-n_d=',np.min(nngb_d),np.median(nngb_d),np.mean(nngb_d),np.std(nngb_d),np.max(nngb_d),' std-log10-n=',np.std(np.log10(nngb_d)))
    print('min/med/mean/std/max-n_g=',np.min(nngb_g),np.median(nngb_g),np.mean(nngb_g),np.std(nngb_g),np.max(nngb_g),' std-log10-n=',np.std(np.log10(nngb_g)))
    print('min/med/mean/std/max-d/g=',np.min(dg),np.median(dg),np.mean(dg),np.std(dg),np.max(dg),' std-log10-n=',np.std(np.log10(dg)))
    print('min/max h_d =',np.min(h_d),np.max(h_d))
    print('min/max h_g =',np.min(h_g),np.max(h_g))
    
    # write the results to an HDF5 binary file
    s0=sdir.split("/"); snapdir_specific=s0[len(s0)-1]; n_s0=1;
    if(len(snapdir_specific)<=1): snapdir_specific=s0[len(s0)-2]; n_s0=2;
    if(snapdir_specific=='output'): snapdir_specific=s0[len(s0)-n_s0-1]; 
    if(len(snapdir_specific)<=1): snapdir_specific=s0[len(s0)-n_s0-2]; 
    ss=snap_ext(snum,four_char=1)
    hmin_string = "%0.0e" % Hmin
    agrain_string = "%0.0d" % a_grain
    #filename='./'+snapdir_specific+'_ngbfile_hmin_'+hmin_string+'_a1e-6_'+agrain_string+'_snap_'+ss+'.dat'
    filename=sdir+'/dust_snap_'+ss+'.hdf5'
    outfi = h5py.File(filename,'w')
    dset_sdir = outfi.create_dataset('Snapshot_Directory',data=sdir)
    dset_snum = outfi.create_dataset('Snapshot_Number',data=snum)
    dset_time = outfi.create_dataset('Snapshot_Time',data=time)
    dset_ndim = outfi.create_dataset('Number_Of_Simulation_Dimensions',data=ndims)
    dset_bsiz = outfi.create_dataset('Periodic_BoxSize',data=BoxSize)
    dset_ndus = outfi.create_dataset('Total_Number_Of_Dust_Particles',data=Ndust)
    dset_ngas = outfi.create_dataset('Total_Number_Of_Gas_Particles',data=Ngas)
    dset_hmin = outfi.create_dataset('Minimum_Smoothing_Length_Allowed',data=Hmin)
    dset_nngb = outfi.create_dataset('Desired_Number_Of_Neighbors',data=Nngb)
    dset___hd = outfi.create_dataset('Smoothing_Length_List_Of_Dust_Neighbors',data=h_d)
    dset___hg = outfi.create_dataset('Smoothing_Length_List_Of_Gas_Neighbors',data=h_g)
    dset___nd = outfi.create_dataset('Mass_Density_List_Of_Dust_Neighbors',data=nngb_d)
    dset___ng = outfi.create_dataset('Mass_Density_List_Of_Gas_Neighbors',data=nngb_g)
    outfi.close()
    # read with:
    #    infi=h5py.File(filename,'r') # (defined same as above)
    #    h_g = np.array(infi["Smoothing_Length_List_Of_Gas_Neighbors"]) # this will act like a normal vector: h_g[:] for all elements
    #    snum = infi["Snapshot_Number"].value # note that for scalars above, we need to add the '.value' to pull out the actual number
    #    infi.close()
    #
    
    return h_d,nngb_d,h_g,nngb_g









# this is the 'master' routine to take a snapshot of a turbulent grain box, pull out the relevant
#   numbers (dust and gas positions), feed it to subroutine get_particle_hsml_ngb which wraps to 
#   the c-libraries that do the neighbor calculations, then saves the outputs to a binary file so 
#   we don't have to re-do the calculation again
#   This version _grid, calclulates density on a Cartesian grid instead. Defaults to using twice as many grid points as gas points, change as desired with resolution_{xyz} keywords
def grain_density_from_snapshot_grid(sdir='/Users/phopkins/Downloads/grains',snum=1,Hmin=0.,boxsize_x=1.,boxsize_y=1.,boxsize_z=1.,DEBUG=False,resolution_x=0,resolution_y=0,resolution_z=0):
    # read in the snapshot data
    Pd=gadget.readsnap(sdir,snum,3); Pg=gadget.readsnap(sdir,snum,0); Ph=gadget.readsnap(sdir,snum,0,header_only=1);
    BoxSize=Ph['boxsize']; time=Ph['time'];
    boxsize_x*=BoxSize; boxsize_y*=BoxSize; boxsize_z*=BoxSize;
    
    
    Ndust = Pd['p'][:,0].size
    Ngas = Pg['p'][:,0].size

    # Setting up grid of search locations
    box_vol=boxsize_x*boxsize_y*boxsize_z;
    res_default=np.int(2*(Ngas/box_vol)**(0.333333))
    twoDQ=False
#    zcoords = np.array()
    if ( np.all(Pg['p'][:,2]==0.)):
        twoDQ=True
        res_default=np.int(2*(Ngas/box_vol)**(0.5))
    
    # JS: Changed to center searching on grid
    if (resolution_x<1):
        resolution_x=np.int(res_default*boxsize_x)
    if (resolution_y<1):
        resolution_y=np.int(res_default*boxsize_y)
    if (resolution_z<1):
        resolution_z=np.int(res_default*boxsize_z)
    
    xgs=np.linspace(0,boxsize_x,num=resolution_x,endpoint=False)+boxsize_x/(resolution_x+1)/2
    ygs=np.linspace(0,boxsize_y,num=resolution_y,endpoint=False)+boxsize_y/(resolution_y+1)/2
    zgs=np.linspace(0,boxsize_z,num=resolution_z,endpoint=False)+boxsize_z/(resolution_z+1)/2
    if twoDQ:
        zgs=0.
    # Grid itself
    if twoDQ:
        xs, ys = np.meshgrid(xgs,ygs,sparse=False, indexing='ij')
        xs = xs.flatten()
        ys = ys.flatten()
        zs = np.zeros(xs.size)
    else:
        xs, ys, zs = np.meshgrid(xgs,ygs,zgs,sparse=False, indexing='ij')
        xs = xs.flatten()
        ys = ys.flatten()
        zs = zs.flatten()

    print('')
    print('Depositing on grid of size',resolution_x,resolution_y,(not twoDQ)*resolution_z,xs.size)
    print('')

    if(DEBUG==True):
        b0=0.2
        ok_d=np.where((Pd['p'][:,0]<b0)&(Pd['p'][:,1]<b0)&(Pd['p'][:,2]<b0))
        ok_g=np.where((Pg['p'][:,0]<b0)&(Pg['p'][:,1]<b0)&(Pg['p'][:,2]<b0))
    else:
        ok_d=np.where(np.abs(Pd['m']) > -1)
        # ok_d=np.where(np.round(Pd['GrainSize']*1.0e6)==a_grain);
        ok_g=np.where(np.abs(Pg['m']) > -1)
    
    DesNgb=0 # Default
    
    
    
    h_d,nngb_d,ndims,Nngb = get_particle_hsml_ngb.get_particle_hsml_ngb(\
                                    xs,ys,zs,Pd['p'][:,0][ok_d],Pd['p'][:,1][ok_d],Pd['p'][:,2][ok_d],\
                                                                        Hmin=Hmin,boxsize_x=boxsize_x,boxsize_y=boxsize_y,boxsize_z=boxsize_z,DesNgb=DesNgb,weight_target=Pd['m'][ok_d])
    h_g,nngb_g,ndims,Nngb = get_particle_hsml_ngb.get_particle_hsml_ngb(\
                                    xs,ys,zs,Pg['p'][:,0][ok_g],Pg['p'][:,1][ok_g],Pg['p'][:,2][ok_g],\
                                                                        Hmin=Hmin,boxsize_x=boxsize_x,boxsize_y=boxsize_y,boxsize_z=boxsize_z,DesNgb=DesNgb,weight_target=Pg['m'][ok_g])
    dg = nngb_d/nngb_g
    print('min/med/mean/std/max-n_d=',np.min(nngb_d),np.median(nngb_d),np.mean(nngb_d),np.std(nngb_d),np.max(nngb_d),' std-log10-n=',np.std(np.log10(nngb_d)))
    print('min/med/mean/std/max-n_g=',np.min(nngb_g),np.median(nngb_g),np.mean(nngb_g),np.std(nngb_g),np.max(nngb_g),' std-log10-n=',np.std(np.log10(nngb_g)))
    print('min/med/mean/std/max-d/g=',np.min(dg),np.median(dg),np.mean(dg),np.std(dg),np.max(dg),' std-log10-n=',np.std(np.log10(dg)))
    print('min/max h_d =',np.min(h_d),np.max(h_d))
    print('min/max h_g =',np.min(h_g),np.max(h_g))
    
    # write the results to an HDF5 binary file
    s0=sdir.split("/"); snapdir_specific=s0[len(s0)-1]; n_s0=1;
    if(len(snapdir_specific)<=1): snapdir_specific=s0[len(s0)-2]; n_s0=2;
    if(snapdir_specific=='output'): snapdir_specific=s0[len(s0)-n_s0-1]; 
    if(len(snapdir_specific)<=1): snapdir_specific=s0[len(s0)-n_s0-2]; 
    ss=snap_ext(snum,four_char=1)
    hmin_string = "%0.0e" % Hmin
    #filename='./'+snapdir_specific+'_ngbfile_hmin_'+hmin_string+'_a1e-6_'+agrain_string+'_snap_'+ss+'.dat'

    # Reshape arrays to cartesian block
    if twoDQ:
        h_d = np.reshape(h_d,(resolution_x,resolution_y))
        h_g = np.reshape(h_g,(resolution_x,resolution_y))
        nngb_d = np.reshape(nngb_d,(resolution_x,resolution_y))
        nngb_g = np.reshape(nngb_g,(resolution_x,resolution_y))
    else:
        h_d = np.reshape(h_d,(resolution_x,resolution_y,resolution_z))
        h_g = np.reshape(h_g,(resolution_x,resolution_y,resolution_z))
        nngb_d = np.reshape(nngb_d,(resolution_x,resolution_y,resolution_z))
        nngb_g = np.reshape(nngb_g,(resolution_x,resolution_y,resolution_z))

    filename=sdir+'/dust_snap_'+ss+'.h5'
    outfi = h5py.File(filename,'w')
    dset_sdir = outfi.create_dataset('Snapshot_Directory',data=sdir)
    dset_snum = outfi.create_dataset('Snapshot_Number',data=snum)
    dset_time = outfi.create_dataset('Snapshot_Time',data=time)
    dset_ndim = outfi.create_dataset('Number_Of_Simulation_Dimensions',data=ndims)
    dset_bsiz = outfi.create_dataset('Periodic_BoxSize',data=BoxSize)
    dset_ndus = outfi.create_dataset('Total_Number_Of_Dust_Particles',data=Ndust)
    dset_ngas = outfi.create_dataset('Total_Number_Of_Gas_Particles',data=Ngas)
    dset_hmin = outfi.create_dataset('Minimum_Smoothing_Length_Allowed',data=Hmin)
    dset_nngb = outfi.create_dataset('Desired_Number_Of_Neighbors',data=Nngb)
    dset___hd = outfi.create_dataset('Smoothing_Length_List_Of_Dust_Neighbors',data=h_d)
    dset___hg = outfi.create_dataset('Smoothing_Length_List_Of_Gas_Neighbors',data=h_g)
    dset___nd = outfi.create_dataset('Mass_Density_List_Of_Dust_Neighbors',data=nngb_d)
    dset___ng = outfi.create_dataset('Mass_Density_List_Of_Gas_Neighbors',data=nngb_g)
    dset___ng = outfi.create_dataset('x_grid',data=xgs)
    dset___ng = outfi.create_dataset('y_grid',data=ygs)
    dset___ng = outfi.create_dataset('z_grid',data=zgs)
    outfi.close()
    # read with:
    #    infi=h5py.File(filename,'r') # (defined same as above)
    #    h_g = np.array(infi["Smoothing_Length_List_Of_Gas_Neighbors"]) # this will act like a normal vector: h_g[:] for all elements
    #    snum = infi["Snapshot_Number"].value # note that for scalars above, we need to add the '.value' to pull out the actual number
    #    infi.close()
    #
    
    return h_d,nngb_d,h_g,nngb_g



def pairCorrelationFunction_2D(x, y, S, rMax, dr):
    """Compute the two-dimensional pair correlation function, also known
    as the radial distribution function, for a set of circular particles
    contained in a square region of a plane.  This simple function finds
    reference particles such that a circle of radius rMax drawn around the
    particle will fit entirely within the square, eliminating the need to
    compensate for edge effects.  If no such particles exist, an error is
    returned. Try a smaller rMax...or write some code to handle edge effects! ;)
    Arguments:
        x               an array of x positions of centers of particles
        y               an array of y positions of centers of particles
        S               length of each side of the square region of the plane
        rMax            outer diameter of largest annulus
        dr              increment for increasing radius of annulus
    Returns a tuple: (g, radii, interior_indices)
        g(r)            a numpy array containing the correlation function g(r)
        radii           a numpy array containing the radii of the
                        annuli used to compute g(r)
        reference_indices   indices of reference particles
    """
    from numpy import zeros, sqrt, where, pi, mean, arange, histogram
    # Number of particles in ring/area of ring/number of reference particles/number density
    # area of ring = pi*(r_outer**2 - r_inner**2)

    # Find particles which are close enough to the box center that a circle of radius
    # rMax will not cross any edge of the box
    bools1 = x > rMax
    bools2 = x < (S - rMax)
    bools3 = y > rMax
    bools4 = y < (S - rMax)
    interior_indices, = where(bools1 * bools2 * bools3 * bools4)
    num_interior_particles = len(interior_indices)

    if num_interior_particles < 1:
        raise  RuntimeError ("No particles found for which a circle of radius rMax\
                will lie entirely within a square of side length S.  Decrease rMax\
                or increase the size of the square.")

    edges = arange(0., rMax + 1.1 * dr, dr)
    num_increments = len(edges) - 1
    g = zeros([num_interior_particles, num_increments])
    radii = zeros(num_increments)
    numberDensity = len(x) / S**2

    # Compute pairwise correlation for each interior particle
    for p in range(num_interior_particles):
        index = interior_indices[p]
        d = sqrt((x[index] - x)**2 + (y[index] - y)**2)
        d[index] = 2 * rMax

        (result, bins) = histogram(d, bins=edges, normed=False)
        g[p, :] = result/numberDensity

    # Average g(r) for all interior particles and compute radii
    g_average = zeros(num_increments)
    for i in range(num_increments):
        radii[i] = (edges[i] + edges[i+1]) / 2.
        rOuter = edges[i + 1]
        rInner = edges[i]
        g_average[i] = mean(g[:, i]) / (pi * (rOuter**2 - rInner**2))

    return (g_average, radii, interior_indices)
####

def pairCorrelationFunction_3D(x, y, z, S, rMax, dr):
    """Compute the three-dimensional pair correlation function for a set of
    spherical particles contained in a cube with side length S.  This simple
    function finds reference particles such that a sphere of radius rMax drawn
    around the particle will fit entirely within the cube, eliminating the need
    to compensate for edge effects.  If no such particles exist, an error is
    returned.  Try a smaller rMax...or write some code to handle edge effects! ;)
    Arguments:
        x               an array of x positions of centers of particles
        y               an array of y positions of centers of particles
        z               an array of z positions of centers of particles
        S               length of each side of the cube in space
        rMax            outer diameter of largest spherical shell
        dr              increment for increasing radius of spherical shell
    Returns a tuple: (g, radii, interior_indices)
        g(r)            a numpy array containing the correlation function g(r)
        radii           a numpy array containing the radii of the
                        spherical shells used to compute g(r)
        reference_indices   indices of reference particles
    """
    from numpy import zeros, sqrt, where, pi, mean, arange, histogram

    # Find particles which are close enough to the cube center that a sphere of radius
    # rMax will not cross any face of the cube
    bools1 = x > rMax
    bools2 = x < (S - rMax)
    bools3 = y > rMax
    bools4 = y < (S - rMax)
    bools5 = z > rMax
    bools6 = z < (S - rMax)

    interior_indices, = where(bools1 * bools2 * bools3 * bools4 * bools5 * bools6)
    num_interior_particles = len(interior_indices)

    if num_interior_particles < 1:
        raise  RuntimeError ("No particles found for which a sphere of radius rMax\
                will lie entirely within a cube of side length S.  Decrease rMax\
                or increase the size of the cube.")

    edges = arange(0., rMax + 1.1 * dr, dr)
    num_increments = len(edges) - 1
    g = zeros([num_interior_particles, num_increments])
    radii = zeros(num_increments)
    numberDensity = len(x) / S**3

    # Compute pairwise correlation for each interior particle
    for p in range(num_interior_particles):
        index = interior_indices[p]
        d = sqrt((x[index] - x)**2 + (y[index] - y)**2 + (z[index] - z)**2)
        d[index] = 2 * rMax

        (result, bins) = histogram(d, bins=edges, normed=False)
        g[p,:] = result / numberDensity

    # Average g(r) for all interior particles and compute radii
    g_average = zeros(num_increments)
    for i in range(num_increments):
        radii[i] = (edges[i] + edges[i+1]) / 2.
        rOuter = edges[i + 1]
        rInner = edges[i]
        g_average[i] = mean(g[:, i]) / (4.0 / 3.0 * pi * (rOuter**3 - rInner**3))

    return (g_average, radii, interior_indices)
    # Number of particles in shell/total number of particles/volume of shell/number density
    # shell volume = 4/3*pi(r_outer**3-r_inner**3)
####
