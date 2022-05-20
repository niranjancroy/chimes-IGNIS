import numpy as np
import h5py as h5py
import os.path
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pylab as pylab
import matplotlib.pyplot as plot
import scipy.interpolate as interpolate
from mpl_toolkits.mplot3d import Axes3D

#
# contains several routines of interest:
#   plotpts_w_gas makes images of gas+dust faceon+edgeon
#     (plotpts_w_gas_loop does this on a loop, for movie-making)
#   plotpts just shows the points
#   gas_rho_image shows gas in color-map, can specify different quantities to map
#   load_snap and other routines check for snapshots and load the hdf5 files
#
#   (note in all of them 'snum' = snapshot number, e.g. 0, 19, 105, etc., 
#     'sdir' = snapshot directory path, can be global or local, 
#     vmin,vmax = min/max values plotted, can be set by hand or if 0, code guesses)
#     'gas_val_toplot' = 'rho' (density), 'b' (b-field strength), 'v' (v-magnitude),
#         'p' (pressure), 'bvec' (b-field with vector lines, overlaid on density in colors), 
#         'vvec' (v-field with vector lines, overlaid on density in colors)
#


def crunchpdfs():
    sdirs=['HII_inner_L0pt1pc/k0_N128','HII_inner_Lau/c0_N128','HII_inner_Lsubau/t0_N128_fixed',
    'HII_outer_Lpc/e0_N128','HII_outer_L200au/s0_N128','HII_outer_Lau/l0_N128',
    'AGB/x0_N128','AGB/u0_N128','AGB/v0_N128','AGB/w0_N128_fixed',
    'CGM_QSO_Lkpc/j0_N128_gamma53_pelecscalecharge_boris','CGM_QSO_Lkpc/j0_N128_boris',
    'HII_inner_L0pt1pc/k0_N128_pelecscalecharge','HII_inner_L0pt1pc/k0_N128_pelecscalecharge_2',
    'jono_default_params/b0_N128','jono_default_params/b0_N128_beta1traditionalbeta','jono_default_params/b0_N64_nocharge','jono_default_params/b0_N128_gamma53',
    'Corona_LRsun/h0_N128_boris','Corona_L2000km/r0_N128_boris','active_runs/m0_N128',
    'WIM_Lpc/f0_N128_boris','WIM_Lpc/f0_N128_5x',
    'WIM_L500au/q0_N128_boris','active_runs/q0_N128_5x','active_runs/g0_N128','active_runs/g0_N128_5x',
    'HII_inner_L0pt1pc/k0_N128_angle','HII_inner_L0pt1pc/k0_N128_angle2',
    'jono_default_params/b0_N64','jono_default_params/b0_N32','jono_default_params/b0_N256_boris',
    'HII_inner_L0pt1pc/k0_N128_mu','HII_inner_L0pt1pc/k0_N128_mu2',
    'jono_default_params/b0_N128_mu','jono_default_params/b0_N128_mu2']
    snames=['k0','c0','t0',
    'e0','s0','l0',
    'x0','u0','v0','w0',
    'j0g53pe','j0',
    'k0pe','k0pe2',
    'b0','b0beta','b0nocharge','b0g53',
    'h0','r0','m0','f0','f05x','q0','q05x','g0','g05x',
    'k0ang1','k0ang2',
    'b0N64','b0N32','b0N256',
    'k0muLO','k0muHI','b0muLO','b0muHI']
    
    sdirs=['active_runs/g0_N128_5x','active_runs/m0_N128','active_runs/g0_N128','active_runs/q0_N128_5x']
    snames=['g05x','m0','g0','q05x']
    
    sdirs=['HII_inner_L0pt1pc/k0_N128_betacorr']; snames=['k0tauHi']
    sdirs=['CGM_QSO_Lkpc/j0_N128_betacorr']; snames=['j0tauHi']
    
    for sdir,sname in zip(sdirs,snames):
        plot_allpdfs(sdir,fname_suffix=sname)
    



def plot_allpdfs(sdir,fname_suffix='',**kwargs):
    pylab.close('all'); plot.figure(1,figsize=(14.,3.25)); f00=14.;
    matplotlib.rcParams.update({'font.size':f00})

    g_wt_vol,g_wt_m,g_v,g_b,g_rho = extractvals(sdir,ptype='PartType0',**kwargs)
    d_wt_vol,d_wt_m,d_v,d_b,d_rho = extractvals(sdir,ptype='PartType3',**kwargs)
    
    nbins=100
    for vplot,wt_v,wt_m,isub,xlab,ylab in zip([g_v,d_v,g_b],[g_wt_vol,d_wt_vol,g_wt_vol],
            [g_wt_m,d_wt_m,g_wt_m],[1,2,3],
            [r'Gas Velocity $({\bf u}_{g}-\langle{\bf u}_{g}\rangle)/c_{s}$',
            r'Dust Velocity $({\bf v}_{d}-\langle{\bf v}_{d}\rangle)/c_{s}$',
            r'Magnetic Field $({\bf B}-\langle{\bf B}\rangle)/\sqrt{4\pi\,P_{0}}$'],
            [r'Fraction','','']):
        pylab.subplot(1,4,isub); pylab.xscale('linear'); pylab.yscale('log'); 
        pylab.xlabel(xlab,fontsize=f00); pylab.ylabel(ylab,fontsize=f00);
        wt_v/=np.sum(wt_v); wt_m/=np.sum(wt_m);
        for k,col,label in zip([0,1,2],['blue','red','green'],[r'$x$',r'$y$',r'$z$']):
            u = 1.*vplot[:,k]; u -= np.sum(u*wt_v);
            y,xb = np.histogram(u,bins=nbins,weights=wt_v,density=True)
            if(isub==1):
                pylab.plot(xb[0:-1]+0.5*np.diff(xb),y/np.max(y),color=col,linewidth=2.,linestyle='-',label=label)
            else:
                pylab.plot(xb[0:-1]+0.5*np.diff(xb),y/np.max(y),color=col,linewidth=2.,linestyle='-')
            y,xb = np.histogram(u,bins=nbins,weights=wt_m,density=True)
            pylab.plot(xb[0:-1]+0.5*np.diff(xb),y/np.max(y),color=col,linewidth=2.,linestyle=':')

            x0=np.linspace(np.min(u),np.max(u),100.); s0=np.sqrt(std_wt(u,wt_v));
            y0=1./np.sqrt(s0*s0*2.*np.pi)*np.exp(-x0*x0/(2.*s0*s0)); ok=(y0>1.e-6)&(y0>np.min(y)); 
            #pylab.plot(x0[ok],y0[ok]/np.max(y0),color=col,linewidth=1.,linestyle='--')
            x00=s0/np.sqrt(2.); y0=np.exp(-np.abs(x0)/x00)/(2.*x00); ok=(y0>1.e-6)&(y0>np.min(y)); 
            #pylab.plot(x0[ok],y0[ok]/np.max(y0),color=col,linewidth=1.,linestyle='-.')
        if(isub==1):
            pylab.legend(loc='upper right',fontsize=f00,handletextpad=0.5,borderpad=0.2,
                frameon=False,columnspacing=0.0,labelspacing=0.0,numpoints=1)
        if(isub==2):
            pylab.plot([0.3e-6,0.3e-6],[-1.e-8,1.e-8],color='black',linewidth=2.,linestyle='-',label=r'$P_{V}$')
            pylab.plot([0.3e-6,0.3e-6],[-1.e-8,1.e-8],color='black',linewidth=2.,linestyle=':',label=r'$P_{M}$')
            pylab.legend(loc='upper right',fontsize=f00,handletextpad=0.5,borderpad=0.2,
                frameon=False,columnspacing=0.0,labelspacing=0.0,numpoints=1)
        pylab.ylim(1.e-6,2.)
        frame1=plot.gca() #get coordinate axes, to use for tick manipulation below
        if(isub>1): frame1.axes.yaxis.set_ticklabels([])
        frame1.axes.tick_params(direction='in',which='both',axis='both',bottom=True,top=True,left=True,right=True)
        frame1.axes.minorticks_on(); 
        frame1.axes.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10.,numticks=55))

    pylab.subplot(1,4,4); pylab.xscale('linear'); pylab.yscale('log'); 
    pylab.xlabel(r'Density $\ln[\rho/\langle\rho\rangle]$',fontsize=f00); pylab.ylabel('',fontsize=f00); 
    for vplot,wt_v,wt_m,col,label in zip([g_rho,d_rho],[g_wt_vol,d_wt_vol],[g_wt_m,d_wt_m],
        ['orange','purple'],[r'Gas',r'Dust']):
            u = 1.*vplot; u /= np.sum(u*wt_v); u = np.log(u);
            print(sdir)
            print('CLUMPING_FACTOR (',label,')=',np.sum(np.exp(u)*wt_m),np.sum(np.exp(2.*u)*wt_v))
            y,xb = np.histogram(u,bins=nbins,weights=wt_v,density=True)
            pylab.plot(xb[0:-1]+0.5*np.diff(xb),y/np.max(y),color=col,linewidth=2.,linestyle='-',label=label)
            y,xb = np.histogram(u,bins=nbins,weights=wt_m,density=True)
            pylab.plot(xb[0:-1]+0.5*np.diff(xb),y/np.max(y),color=col,linewidth=2.,linestyle=':')

            x0=np.linspace(np.min(u),np.max(u),100.); s0=np.sqrt(std_wt(u,wt_v)); dx=x0+s0*s0/2.
            y0=1./np.sqrt(s0*s0*2.*np.pi)*np.exp(-dx*dx/(2.*s0*s0)); ok=(y0>1.e-6)&(y0>np.min(y)); 
            #pylab.plot(x0[ok],y0[ok]/np.max(y0),color=col,linewidth=1.,linestyle='--')
            x00=s0/np.sqrt(2.); dx=x0 - np.log(1.-x00*x00); 
            y0=np.exp(-np.abs(dx)/x00)/(2.*x00); ok=(y0>1.e-6)&(y0>np.min(y)); 
            #pylab.plot(x0[ok],y0[ok]/np.max(y0),color=col,linewidth=1.,linestyle='-.')
    pylab.legend(loc='lower left',fontsize=f00,handletextpad=0.5,borderpad=0.0,
        frameon=False,columnspacing=0.0,labelspacing=0.0,numpoints=1)
    pylab.ylim(1.e-6,2.)
    frame1=plot.gca() #get coordinate axes, to use for tick manipulation below
    frame1.axes.yaxis.set_ticklabels([])
    frame1.axes.tick_params(direction='in',which='both',axis='both',bottom=True,top=True,left=True,right=True)
    frame1.axes.minorticks_on(); 
    frame1.axes.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10.,numticks=55))

    pylab.subplots_adjust(left=0.1,bottom=0.1,right=1-0.01,top=1-0.01,hspace=0.01,wspace=0.01);
    pylab.savefig('pdf_plot_'+fname_suffix+'.pdf',transparent=True,bbox_inches='tight',pad_inches=0)


def extractvals(sdir, ptype='PartType0',
    sdir_sub='/output/',
    sdir_master='/panfs/ds08/hopkins/phopkins/data2/dust_instability_boxes/mhd_rhi_dust_boxes/',
    n_touse=10):
    
    sdir=sdir_master+'/'+sdir+'/'+sdir_sub
    snums=np.array([],dtype='int')
    for snum in range(800):
        fname,fname_base,fname_ext = check_if_filename_exists(sdir,snum)
        if(fname != 'NULL'): snums=np.append(snums,snum)
    if(snums.size < 10*n_touse): n_touse=snums.size/10
    if(n_touse < 1): n_touse=1;
    snums=snums[snums.size-n_touse::]
    dv_g_vwt=np.zeros(3); dv_g_mwt=np.zeros(3); dv_d_vwt=np.zeros(3); dv_d_mwt=np.zeros(3); dB_g_vwt=np.zeros(3); dB_g_mwt=np.zeros(3);
    drho_g_vwt=0.; drho_g_mwt=0.; drho_d_vwt=0.; drho_d_mwt=0.; drho_dg_vwt=0.; drho_dg_mwt=0.;
    print(snums)
    a_wtv=np.zeros(0); a_wtm=np.zeros(0); a_v=np.zeros((0,3)); a_rho=np.zeros(0); a_b=np.zeros((0,3))
    for snum in snums:
        print('scanning snapshot ',snum)
        F=load_snap(sdir,snum);
        print('Time == ',F['Header'].attrs['Time'])
        P=F[ptype]; v=np.array(P['Velocities']); m=np.array(P['Masses']);
        if(ptype=='PartType0'):
            rho=np.array(P['Density']); b=np.array(P['MagneticField']);
        else:
            rho,h_d,rho_gn,h_gn=load_snap_dustrho(sdir,snum); b=np.zeros((0,3));
        rho+=1.e-10; m[(m<=0)] += m[(m>0.)].min(); rho[np.isfinite(rho)==False]=rho[np.isfinite(rho)==True].min()
        vol=m/rho; vol/=np.sum(vol); m/=np.sum(m);
        a_wtv = np.concatenate((a_wtv,vol))
        a_wtm = np.concatenate((a_wtm,m))
        v0=np.zeros(3); b0=np.zeros(3);
        for k in [0,1,2]:
            v0[k]=np.sum(v[:,k]*vol)
            v[:,k]-=v0[k];
            if(ptype=='PartType0'):
                b0[k]=np.sum(b[:,k]*vol);
                b[:,k]-=b0[k];
        a_v = np.concatenate((a_v,v))
        a_b = np.concatenate((a_b,b))
        rho/=np.sum(rho*vol)
        a_rho = np.concatenate((a_rho,rho))
        F.close()
    return a_wtv,a_wtm,a_v,a_b,a_rho





def plotpts_w_gas(snum=0,sdir='./output',ptype='PartType3',width=0.05,cut_dust=1.,alpha=0.1,markersize=5.,
        vmin=0,vmax=0,forsavedfigure=False,gas_val_toplot='rho',ptype_im='PartType0',
        zmed_set=-1.e10,cmap='terrain'):
    pylab.close('all'); 
    plot.figure(1,figsize=(21.,7.))
    P_File=load_snap(sdir,snum);
    print('Time == ',P_File['Header'].attrs['Time'])
    P = P_File[ptype]
    print('Dust-to-Gas Mass Ratio = ',np.sum(P_File['PartType3']['Masses'][:])/np.sum(P_File['PartType0']['Masses'][:]))
    a = P['GrainSize'][:]
    print('Grain Size: min=',a.min(),' max=',a.max())
    Pc=np.array(P['Coordinates']);
    vx = np.array(P['Velocities'][:,0])
    vy = np.array(P['Velocities'][:,1])
    vz = np.array(P['Velocities'][:,2])
    print('Var v_x_dust == ',np.min(vx),np.median(vx),np.max(vx),np.std(vx),' (min/median/max/std)')
    print('Var v_y_dust == ',np.min(vy),np.median(vy),np.max(vy),np.std(vy),' (min/median/max/std)')
    print('Var v_z_dust == ',np.min(vz),np.median(vz),np.max(vz),np.std(vz),' (min/median/max/std)')
    
    for subplot,xz,yz in zip([1,2,3],[0,1,0],[0,0,1]): 
        xplotc=0; yplotc=1; depth_c=2; #plot xy, z=depth
        if(xz==1): 
            xplotc=0; yplotc=2; depth_c=1; #plot xz, y=depth
        if(yz==1): 
            xplotc=1; yplotc=2; depth_c=0; #plot yz, x=depth
        
        quiet=True; 
        if(subplot==1): quiet=False # give numbers once
        z=Pc[:,depth_c]; zmx=np.max(z)-np.min(z); z0=np.median(z); 
        if(zmed_set > -1.e9): z0=zmed_set 
        zz=np.abs(z-z0); zz[(zz>0.5*zmx)] = zmx-zz[(zz>0.5*zmx)]
        ok=np.where(zz < width)
        x=Pc[:,xplotc][ok]; y=Pc[:,yplotc][ok]; 
        pylab.subplot(1,3,subplot)
        ok_r = (np.random.rand(x.size) < cut_dust)
        gas_rho_image(snum=snum,sdir=sdir,xmax=1.,xz=xz,yz=yz,gas_val_toplot=gas_val_toplot,zmed_set=z0,
            vmin=vmin,vmax=vmax,quiet=quiet,cmap=cmap,ptype=ptype_im)
        if(forsavedfigure==True):
            pylab.plot(x[ok_r],y[ok_r],marker='.',markersize=markersize,alpha=alpha,
                linestyle='',color='black',markeredgecolor='None',zorder=3)
        else:
            pylab.plot(x[ok_r],y[ok_r],marker='.',markersize=markersize,alpha=alpha,
                linestyle='',color='black',markeredgecolor='None',zorder=3)
        frame1=plot.gca() #get coordinate axes, to use for tick manipulation below
        frame1.axes.xaxis.set_ticklabels([]) #no tick names
        frame1.axes.yaxis.set_ticklabels([]) #no tick names
        frame1.axes.get_xaxis().set_ticks([]) #no ticks
        frame1.axes.get_yaxis().set_ticks([]) #no ticks

    if(forsavedfigure==True): 
        ext='000'+str(snum);
        if (snum>=10): ext='00'+str(snum)
        if (snum>=100): ext='0'+str(snum)
        if (snum>=1000): ext=str(snum)
        pylab.savefig('im_'+ext+'.png',dpi=150,bbox_inches='tight',pad_inches=0)
    P_File.close()



def symbology(sdir):
    marker='^'; color='black'; linestyle='-.'; markersize=6.; linewidth=1.5; 
    if('AGB' in sdir or 'w0_' in sdir or 'u0_' in sdir or 'v0_' in sdir or 'x0_' in sdir):
        color='blue'
    if('HII_inner' in sdir or 't0_' in sdir or 'c0_' in sdir or 'k0_' in sdir or 'HII-near' in sdir):
        color='orange'
    if('HII_outer' in sdir or 'l0_' in sdir or 's0_' in sdir or 'e0_' in sdir or 'HII-far' in sdir):
        color='green'
    if('WIM' in sdir or 'g0_' in sdir or 'q0_' in sdir or 'f0_' in sdir):
        color='violet'
    if('Corona' in sdir or 'm0_' in sdir or 'r0_' in sdir or 'h0_' in sdir):
        color='brown'
    if('CGM' in sdir or 'j0_' in sdir):
        color='red'
    if('jono' in sdir or 'b0_' in sdir or 'Example' in sdir):
        color='black'

    # 'S' boxes
    if('w0_' in sdir or 't0_' in sdir or 'l0_' in sdir or 'g0_' in sdir or 'm0_' in sdir):
        marker='o'; linestyle=':'; linewidth=2.5; 
    # 'M' boxes
    if('j0_' in sdir or 'b0_' in sdir or 'v0_' in sdir or 'c0_' in sdir or 's0_' in sdir or 'q0_' in sdir or 'r0_' in sdir):
        marker='s'; linestyle='-'; linewidth=1.5; 
    # 'L' boxes
    if('u0_' in sdir or 'k0_' in sdir or 'e0_' in sdir or 'f0_' in sdir or 'h0_' in sdir):
        marker='*'; linestyle='--'; linewidth=2.0; markersize*=1.5
    # special cases
    if('tau' in sdir or '1x' in sdir or 'beta1' in sdir or 'gamma' in sdir or 'pe' in sdir or 'mu' in sdir or 'ang' in sdir):
        markersize*=0.5
    
    return marker,markersize,color,linestyle,linewidth
          



def gas_veldisp_vs_time(dodust=False):
    pylab.close('all'); plot.figure(1,figsize=(6.,4.5)); f00=16.
    matplotlib.rcParams.update({'font.size':f00})

    ptype='PartType0'; pylab.axis([0.,14.95,1.1e-4,50.]); fname='gasvel_vs_time'; ylab=r'Velocity Dispersion: $|\delta {\bf u}_{g}|$ / $c_{s}^{0}$'
    if(dodust==True):
        ptype='PartType3'; pylab.axis([0.,14.95,1.1e-4,50.]); fname='dustvel_vs_time'; ylab=r'Velocity Dispersion: $|\delta {\bf v}_{d}|$ / $c_{s}^{0}$'

    xlab=r'Time / $t_{s}^{0}$'
    xlab=r'Time: t / $\langle t_{\rm grow}^{0}[L_{\rm box}] \rangle$'
    sdirs=['jono_default_params/b0_N128',
    'AGB/x0_N128','AGB/u0_N128','AGB/v0_N128','AGB/w0_N128_fixed',
    'HII_inner_L0pt1pc/k0_N128','HII_inner_Lau/c0_N128','HII_inner_Lsubau/t0_N128_fixed',
    'HII_outer_Lpc/e0_N128','HII_outer_L200au/s0_N128','HII_outer_Lau/l0_N128',
    'WIM_Lpc/f0_N128_5x','active_runs/q0_N128_5x','active_runs/g0_N128_5x',
    'Corona_LRsun/h0_N128_boris','Corona_L2000km/r0_N128_boris','active_runs/m0_N128',
    'CGM_QSO_Lkpc/j0_N128_boris']
    smaster='/panfs/ds08/hopkins/phopkins/data2/dust_instability_boxes/mhd_rhi_dust_boxes/'
    tmin =[1.5, 
           2.e-4,           4.e-3,         5.e-2,        1.0, 
           1.e-2,           0.7,            10.,    
           4.e-3,           2.e-2,          0.7,       
            1.e-3,          1.4e-3,         3.e-2,
            10.e-2,          1.0,            14.,
           25.]
    tmax =[9.0,
           3.e-3,           3.e-2,         6.e-1,        4.0,
           0.25,            5.0,            50.,    
           1.2e-2,          2.e-1,          3.,         
           1.e-3,           1.4e-2,         1.0,
           0.6,             10.,            140.,
           3000.]
    tbox =[0.34,
            3.e5 ,           920.,           3.1 ,        0.01,
            1000.,          2.6,            0.0088,     
            330.,           1.0,            0.0032,       
            100.,           0.21,           4.8e-4,
            20.,            0.067,          3.3e-4,
            0.29]



    tmin=np.array(tmin); tmax=np.array(tmax); tbox=np.array(tbox);
    omega_eff = np.sqrt(tmin*tmax)
    omega_eff = tmin
    tgrow_code = 1./(tbox * omega_eff);

    snapdirs=np.array(sdirs); snapdir_names=snapdirs;
    colors=['black','red','blue','green','orange','magenta','lime','purple','yellow','black','red','blue','green','orange','magenta','lime','purple','yellow']
    linest=['-',':','--','-.','-','--','-.',':','-','--','-.',':','-','--','-.',':','-',':','--','-.','-','--','-.',':','-','--','-.',':','-','--','-.',':']
    for sdir,label,color,lsty,tnorm in zip(snapdirs,snapdir_names,colors,linest,tgrow_code):    
        sdirX=smaster+sdir+'/output/'
        print('snapdir = ',sdir)
        snums=np.array([],dtype='int')
        for snum in range(800):
            fname_x,fname_base,fname_ext = check_if_filename_exists(sdirX,snum)
            if(fname_x != 'NULL'): snums=np.append(snums,snum)
        v2=np.zeros(snums.size); t=np.zeros(snums.size);
        for snum,i in zip(snums,range(snums.size)):
            print(snum)
            F=load_snap(sdirX,snum); 
            v=np.array(F[ptype+'/Velocities']); m=np.array(F[ptype+'/Masses']); 
            t[i]=F['Header'].attrs['Time'];
            wt=m/np.sum(m);
            vt=v[:,0]; v0=np.sum(vt*wt); v2[i]+=np.sum(vt*vt*wt)-v0*v0;
            vt=v[:,1]; v0=np.sum(vt*wt); v2[i]+=np.sum(vt*vt*wt)-v0*v0;
            vt=v[:,2]; v0=np.sum(vt*wt); v2[i]+=np.sum(vt*vt*wt)-v0*v0;
            F.close()
        v2_min = 1.e-10
        v2[(np.isfinite(v2)==False)]=v2_min
        v2[(v2 <= v2_min)]=v2_min;
        if(tnorm <= 0.): tnorm=np.max(t);
        marker,markersize,color,linestyle,linewidth = symbology(sdir)
        xplot = t/tnorm; yplot = np.sqrt(v2);
        ok=(xplot < 15.5); xplot=xplot[ok]; yplot=yplot[ok];
        pylab.plot(xplot,yplot,color=color,linestyle=linestyle,linewidth=linewidth)
        dx = np.median(np.abs(np.diff(xplot)))
        nsamp = np.round(0.25/dx).astype('int')
        if(nsamp<=1): nsamp=1;
        pylab.plot(xplot[0::nsamp],yplot[0::nsamp],markersize=markersize,
            marker=marker,markerfacecolor='None',color=color,linestyle='')
        print(t, np.sqrt(v2))
    #pylab.legend(loc='lower right',fontsize=f00,handletextpad=0.5,borderpad=0.5,
    #    frameon=False,columnspacing=0.0,labelspacing=0.0);

    pylab.xlabel(xlab)
    pylab.ylabel(ylab)
    pylab.xscale('linear'); pylab.yscale('log');
    frame1=plot.gca() #get coordinate axes, to use for tick manipulation below
    frame1.axes.tick_params(direction='in',which='both',axis='both',
        bottom=True,top=True,left=True,right=True)
    frame1.axes.minorticks_on(); 
    #pylab.xscale('log'); pylab.yscale('log'); 
    xmin=3.; ymin=2.e-4;
    x0=np.linspace(xmin,xmin+2.,10); y0=ymin*np.exp(x0-np.min(x0)); pylab.plot(x0,y0,color='magenta',linestyle=':',linewidth=2.)
    x0=np.linspace(xmin,xmin+2./10.,10); y0=ymin*np.exp(10.*(x0-np.min(x0))); pylab.plot(x0,y0,color='magenta',linestyle=':',linewidth=2.)
    pylab.axis([0.,15.,1.1e-4,50.])
    #pylab.axis([1.e-4,600.,1.1e-4,50.]); pylab.xscale('log')

    pylab.subplots_adjust(left=0.1,bottom=0.1,right=1-0.01,top=1-0.01,hspace=0.01,wspace=0.01);
    pylab.savefig(fname+'.pdf',transparent=True,bbox_inches='tight',pad_inches=0)
        


def plot_corr(activelist=['b0_',
    'b0_beta1','b0_gamma53','b0_mu0pt001','b0_mu0pt1',
    'w0_','v0_','u0_','x0_',
    't0_','c0_','k0_','k0_tauHi','k0_pe','k0_pe2','k0_ang1','k0_ang2','k0_mu0pt001','k0_mu0pt1',
    'l0_','s0_','e0_',
    'm0_','r0_','h0_','h0_tauLO',
    'g0_1x','q0_1x','f0_1x',
    'g0_5x','q0_5x','f0_5x',
    'j0_','j0_pegamma53','j0_tauLO','j0_tauLOpegamma53']):
    pylab.close('all'); plot.figure(1,figsize=(10.,5.5)); f00=12.;
    matplotlib.rcParams.update({'font.size':f00})

    mastervals=[
    ['b0_' , 0.075 , 0.025 , 0.017 , 0.086 , 0.063 , 0.32 , 0.075 , 0.025 , 6.4e-3 , 0.013 , 1.2 , 1.2, 0.89, 25., 29., 50., 0.05, 0.34, 5., 1., 1.5, 8.0, 0.01],
    ['b0_beta1', 0.067 , 0.015 , 6.6e-3 , 0.095 , 0.086 , 0.20 , 0.067 , 0.013 , 3.1e-3 , 3.7e-3 , 1.7 , 1.7, 0.89, 25., 29., 50., 0.05, 0.34, 5., 0.5, 1.5, 8.0 , 0.01],
    ['b0_gamma53' , 0.088 , 0.039 , 0.020 , 0.12 , 0.11 , 0.33 , 0.089 , 0.037 , 0.013 , 8.2e-3 , 1.2, 1.2, 0.89, 25., 29., 50., 0.05, 0.34, 5., 1., 1.5, 8.0, 0.01],
    ['b0_mu0pt001' , 4.3e-3 , 1.8e-3, 6.7e-3 , 2.6e-3, 3.2e-3 , 0.035, 5.1e-3, 1.8e-3 ,4.1e-3, 4.9e-3, 2.7, 2.7, 0.89, 25., 29., 50., 0.05, 0.34, 5., 1., 0.35, 3.1, 0.001],
    ['b0_mu0pt1' ,   0.18, 0.15 , 0.11 , 0.32, 0.30, 0.41, 0.17, 0.12, 0.056, 0.055, 0.89, 0.88, 0.89, 25., 29., 50., 0.05, 0.34, 5., 1., 5.0, 18.0, 0.1],
    ['w0_' , 0.052 , 0.050 , 0.014 , 6.9e-4, 9.4e-4, 7.6e-4 , 0.0086, 0.0094, 0.049 , 0.052 , 0.75 , 0.74, 3., 8.4, 1.e-3, 2.8e-3, 0.707, 0.01, 270., 1. , 1., 4., 0.01],
    ['v0_', 0.031 , 0.010 , 0.038 , 0.10 , 6.8e-3 , 0.11 , 0.066 , 4.2e-3 , 0.071 , 0.074 , 0.90 , 0.88, 3., 8.4, 1.e-3, 2.8e-3, 0.707, 3.1, 0.9, 1. , 5.e-2, 6.e-1, 0.01],
    ['u0_', 1.0 , 0.50 , 1.2 , 1.7 , 0.50 , 1.9 , 1.2 , 0.51 , 1.2 , 0.85 , 0.88 , 0.76, 3., 8.4, 1.e-3, 2.8e-3, 0.707, 920., 0.003, 1. , 4.e-3, 3.e-2, 0.01],
    ['x0_', 19 , 4.3 , 18 , 20 , 3.2 , 20 , 11 , 2.6 , 12 , 1.1 , 1.1 , 0.70, 3., 8.4, 1.e-3, 2.8e-3, 0.707, 3.e5, 1.e-5, 1. , 2.e-4, 3.e-3, 0.01],
    ['t0_', 0.06 , 0.057 , 0.013 , 3.5e-3 , 3.6e-3 , 3.7e-3 , 5.1e-3 , 4.7e-3 , 0.024 , 0.08 , 1.2 , 1.2, 4., 20., 2.1, 24., 0.707, 0.0088, 390., 10., 10., 50., 0.01] ,
    ['c0_', 0.15 , 0.13 , 0.20 , 0.76 , 0.80 , 0.65 , 0.15 , 0.11 , 0.14 , 0.058 , 1.4 , 1.4, 4., 20., 2.1, 24., 0.707, 2.6, 1.3, 10. , 0.7, 5.0, 0.01],
    ['k0_', 3.3 , 1.8 , 3.6 , 13 , 15 , 12 , 3.7 , 1.9 , 3.5 , 1.2 , 1.5 , 1.3 , 4., 20., 2.1, 24., 0.707, 1000., 0.0034, 10., 1.e-2, 0.25, 0.01],
    ['k0_tauHi', 4.7 , 1.6 , 4.8 , 10. , 3.0 , 10. , 4.1 , 1.9 , 4.0 , 1.3 , 2.0 , 1.9 , 4., 20., 7., 24., 0.707, 1000., 0.0034, 10., 1.e-2, 0.25, 0.01],
    ['k0_pe', 2.9 , 2.6 , 2.7 , 29 , 30 , 26 , 2.5 , 1.6 , 2.6, 1.6 , 1.8 , 1.7, 4., 20., 2.1, 24., 0.707, 1000., 0.0034, 10., 1.e-2, 0.25, 0.01],
    ['k0_pe2', 6.9, 7.0, 6.3, 80, 81, 63, 4.3, 4.0, 5.0 , 2.0, 2.2 , 1.4, 4., 20., 2.1, 24., 0.70, 1000., 0.0034, 10., 1.e-2, 0.25, 0.01] ,
    ['k0_ang1', 7.2, 2.1, 4.7, 15, 3.4, 8.5, 4.8, 2.2, 3.0 , 1.2, 2.2 , 2.2, 4., 20., 2.1, 24., 0.70, 1000., 0.0034, 10., 1.5e-2, 0.40, 0.01] ,
    ['k0_ang2', 2.8, 1.5, 7.1, 4.4, 2.7, 10.7, 3.0, 1.7, 6.8, 1.3, 1.9, 1.8, 4., 20., 2.1, 24., 0.70, 1000., 0.0034, 10., 0.5e-2, 0.125, 0.01] ,
    ['k0_mu0pt001', 1.1, 0.36, 1.2, 3.0, 1.2, 3.2, 0.88, 0.28, 0.95, 0.56, 2.4, 2.4, 4., 20., 2.1, 24., 0.70, 1000., 0.0034, 10., 0.3e-2, 0.07, 0.001] ,
    ['k0_mu0pt1', 14., 3.4, 14., 15., 7.4, 15., 12., 3.5, 12.4, 1.7, 1.7, 1.1, 4., 20., 2.1, 24., 0.70, 1000., 0.0034, 10., 3.e-2, 0.9, 0.1] ,
    ['l0_' , 1.8e-3 , 1.8e-3 , 2.0e-3 , 4.7e-5 , 1.4e-5 , 2.1e-5 , 1.4e-3 , 1.4e-3 , 1.4e-3 , 4.3e-4 , 0.63 , 0.62, 0.15, 0.45, 4., 16., 0.707, 0.0032, 500., 10. , 0.7, 7.0, 0.01] ,
    ['s0_', 2.9e-3 , 2.7e-3 , 4.5e-3 , 1.5e-3 , 1.6e-3 , 3.5e-3 , 2.4e-3 , 1.9e-3 , 2.3e-3 , 7.4e-4 , 1.3 , 1.3, 0.15, 0.45, 4., 16., 0.707, 1.0, 1.7, 10. , 2.e-3, 2.e-1, 0.01],
    ['e0_', 0.55 , 0.20 , 0.58 , 1.1 , 0.90 , 1.0 , 0.43 , 0.19 , 0.44 , 0.28 , 1.5 , 1.5, 0.15, 0.45, 4., 16., 0.707, 330., 0.005, 10. , 4.e-3, 1.2e-2, 0.01],
    ['g0_1x', 9.8e-4 , 9.8e-4 , 1.1e-3 , 6.1e-4 , 6.1e-4 , 9.0e-6 , 9.7e-4 , 9.7e-4 , 8.3e-4 , 7.2e-4 , 0.46 , 0.45 , 0.05, 0.1, 100., 160., 0.707, 4.8e-4, 3400., 1., 0.8e-2, 0.7, 0.01],
    ['g0_5x', 9.8e-4 , 9.8e-4 , 1.1e-3 , 6.4e-4 , 6.4e-4 , 9.4e-6 , 9.7e-4 , 9.7e-4 , 1.1e-3 , 7.2e-4 , 0.49 , 0.50 , 0.05, 0.1, 100., 160., 0.707, 4.8e-4, 3400., 1., 3.e-2, 1.0, 0.01],
    ['q0_1x', 10.e-4 , 10.e-4 , 1.1e-3 , 0.6e-3 , 0.6e-3 , 0.2e-3 , 9.7e-4 , 9.7e-4 , 8.5e-4 , 7.3e-4 , 0.25 , 0.22, 0.05, 0.1, 100., 160., 0.707, 0.067, 240., 1. ,3.0e-4,4.0e-3, 0.01],
    ['q0_5x', 9.9e-4 , 9.8e-4 , 1.1e-3 , 1.4e-3 , 1.4e-3 , 1.2e-4 , 9.7e-4 , 9.7e-4 , 9.4e-4 , 7.2e-4 , 0.36 , 0.073, 0.05, 0.1, 100., 160., 0.707, 0.067, 240., 1. ,1.4e-3,1.4e-2, 0.01],
    ['f0_1x', 1.5e-3 , 1.3e-3 , 2.3e-3 , 1.1e-3 , 7.3e-3 , 1.2e-3 , 1.1e-3 , 1.1e-3 , 8.0e-4 , 8.0e-4 , 0.25 , 0.24, 0.05, 0.1, 100., 160., 0.707, 20., 0.8, 1. , 1.e-3, 1.e-3, 0.01],
    ['f0_5x', 0.032  , 0.021  , 0.088  ,   5.30 ,   5.30  , 0.085 ,  3.6e-3, 3.6e-3, 3.9e-3 ,  3.5e-3 , 0.37  , 0.35, 0.05, 0.1, 100., 160., 0.707, 20., 0.8, 1. , 1.e-3, 1.e-3, 0.01],
    ['m0_', 0.011 , 0.011 , 0.013 , 0.012 , 0.012 , 3.5e-5 , 0.011 , 0.011 , 0.013 , 0.060 , 0.63 , 0.59, 20., 480., 3400., 1700., 0.707, 3.3e-4, 4.8e4, 0.001, 14., 140., 0.01],
    ['r0_', 0.012 , 0.012 , 0.014 , 0.62 , 0.63 , 0.01 , 0.011 , 0.011 , 0.022 , 0.060 , 0.86 , 0.86, 20., 480., 3400., 1700., 0.707, 0.067, 240., 0.001 , 1.0, 10., 0.01],
    ['h0_', 0.27 , 0.24 , 0.086 , 24 , 24 , 3.1 , 0.076 , 0.069 , 0.15 , 0.066 , 0.57 , 0.57, 20., 480., 3400., 1700., 0.707, 20., 0.8, 0.001 , 2.e-2, 0.2, 0.01],
    ['h0_tauLO', 0.14 , 0.05 , 0.21 , 0.17 , 0.13 , 0.21 , 0.023, 0.023, 0.013 , 0.085 , 0.60 , 0.57, 20., 480., 100., 1700., 0.707, 20., 0.8, 0.001 , 2.e-2, 0.2, 0.01],
    ['j0_tauLO', 0.031 , 0.031 , 0.036 , 2.0 , 2.0 , 1.9 , 0.024 , 0.021 , 0.019 , 0.021 , 0.27 , 0.26, 9.5, 100., 4300., 1.e6, 0.707, 0.29, 25., 1000. , 25., 3000., 0.01],
    ['j0_tauLOpegamma53', 0.047 , 0.046 , 0.057 , 2.9 , 2.9 , 2.9 , 0.035 , 0.032 , 0.023 , 0.025 , 0.27 , 0.25 , 9.5, 100., 4300., 1.e6, 0.707, 0.29, 25., 1000. , 25., 3000. , 0.01],
    ['j0_', 0.83 , 0.82 , 0.63 , 23. , 23. , 22. , 0.04 , 0.04 , 0.042 , 0.34 , 0.28 , 0.34, 9.5, 100., 1.4e5, 1.e6, 0.707, 0.29, 25., 1000. , 25., 3000., 0.01],
    ['j0_pegamma53', 2.5, 2.5 , 2.0 , 63. , 63., 63. , 0.14, 0.14, 0.13 , 0.56 , 0.25, 0.54 , 9.5, 100., 1.4e5, 1.e6, 0.707, 0.29, 25., 1000. , 25., 3000. , 0.01]
    ]
    
    # first set up the 'symbology' legends
    slist=np.array([r'Example',r'AGB',r'HII-near',r'HII-far',r'WIM',r'Corona',r'CGM'])
    for sdir,i in zip(slist,range(slist.size)):
        pylab.subplot(2,3,1)
        frame1=plot.gca() #get coordinate axes, to use for tick manipulation below
        marker,markersize,color,linestyle,linewidth = symbology(sdir)
        pylab.text(0.065,0.95-0.08*i,sdir,color=color,ha='left',va='top',
            transform=frame1.transAxes,fontsize=f00)
    pylab.subplot(2,3,2)
    names=np.array([r'S',r'M',r'L',r'XL']); slist=np.array(['w0_','v0_','u0_','x0_'])
    for sdir,i,name in zip(slist,range(slist.size),names):
        marker,markersize,color,linestyle,linewidth = symbology(sdir)
        pylab.plot([-1.,-1.],[-1.,-1],markersize=0.5*markersize,label=name,
            marker=marker,markerfacecolor='None',color='black',linestyle='',linewidth=2.)
    pylab.legend(loc='upper left',fontsize=f00,handletextpad=0.,borderpad=0.2,
        frameon=False,columnspacing=0.0,labelspacing=0.0,numpoints=1)


    for key,ipanel in zip(['gasvel','bmag','egy','gasden','dustvel','dustden'],[1,2,3,4,5,6]):
        pylab.subplot(2,3,ipanel)
        x=np.zeros(0); y=np.zeros(0);
        for isim in range(np.shape(mastervals)[0]):
            vals=mastervals[isim]
            sim=vals[0]
            if(sim in activelist):
                marker,markersize,color,linestyle,linewidth = symbology(sim)
                ux=vals[1]; uy=vals[2]; uz=vals[3]; u=np.sqrt(ux*ux+uy*uy+uz*uz); # gas vel
                vx=vals[4]; vy=vals[5]; vz=vals[6]; v=np.sqrt(vx*vx+vy*vy+vz*vz); # dust vel
                bx=vals[7]; by=vals[8]; bz=vals[9]; b=np.sqrt(bx*bx+by*by+bz*bz); # B
                rhog=vals[10]; rhod=vals[11]; rho_dg=vals[12]; # rho_gas, rho_dust, rho_dust/rho_gas
                ws=vals[13]; a_tilde=vals[14]; tau=vals[15]; phi_tilde=vals[16]; cos_ba=vals[17]; Lbox_csts=vals[18]; alpha_tilde=vals[19]; beta=vals[20]; # run initial conditions
                omega_ts_kbox=vals[21]; omega_ts_kmax=vals[22]; mu=vals[23]; # analytic maximum growth rate at box scale and resolution scale, respectively; dust-to-gas ratio
                vA_0=np.sqrt(1./beta); B_final=np.sqrt(1./beta + b*b); vA=B_final; vfast=np.sqrt(1.+vA*vA);
                vA_perp=b; vfast_perp=np.sqrt(1.+vA_perp*vA_perp);

                # for S [circles]; omega_L ok-ish. B_final bad, "b" better but bad for all others; dif't omega-weight ok, but not better than default omega-weight; 
                # model for egy or 'u' matching sqrt[a] is surprisingly good

                if('gasvel' in key):
                    yplot=u; xplot=omega_ts_kbox * Lbox_csts;
                    #xplot=np.sqrt(omega_ts_kmax*omega_ts_kbox)*Lbox_csts / 5. # sort-of-weighted average of (omega[k]/k), which is what matters
                    #yplot=u; xplot=np.sqrt(0.01*a_tilde/alpha_tilde);
                    # CGM outlier here looks good if plot vdust vs this, makes sense?

                    yplot=u; xplot=omega_ts_kbox * Lbox_csts/2.;

                if('bmag' in key):
                    yplot=b; xplot=mu*a_tilde/alpha_tilde/B_final; # agreement here basically defines where fields bent vs fields coherent

                    #yplot=v; xplot=omega_ts_kbox * Lbox_csts / 2.;
                    #xplot=np.sqrt(omega_ts_kmax*omega_ts_kbox)*Lbox_csts / 5. # sort-of-weighted average of (omega[k]/k), which is what matters


                if('egy' in key):
                    xplot=u*u/2.; yplot=b*b/2.; 
                    
                    #yplot=(bx*bx + by*by + B_final*bz)/2.
                    
                    #yplot=u;
                    #xplot=np.sqrt(mu*a_tilde/alpha_tilde * np.sqrt(omega_ts_kmax*omega_ts_kbox)*Lbox_csts / 5.)
                    #xplot=np.sqrt(mu*a_tilde/alpha_tilde * np.sqrt(omega_ts_kmax*omega_ts_kbox)*Lbox_csts ) / 0.1

                    #xplot=(mu*a_tilde/alpha_tilde*Lbox_csts * np.sqrt(omega_ts_kmax*omega_ts_kbox) * ws)**(1./3.)

                    #omega = np.sqrt(omega_ts_kmax*omega_ts_kbox)
                    #cut = ws * np.sqrt(0.01/omega)
                    #if(u > cut):
                    #    print 'u > cut',sim,u,cut
                    #    xplot=np.sqrt(omega_ts_kmax*omega_ts_kbox)*Lbox_csts / 5.
                    #else:
                    #    print 'u < cut',sim,u,cut
                    #    xplot=(0.01*Lbox_csts*ws*ws)**(1./3.)

                if('gasden' in key):
                    yplot=rhog; xplot=u;
                    yplot=rhog; xplot=np.sqrt(np.log(1.+(u/3.)**2.));
                    
                    #yplot=u
                    #xplot=(0.01*Lbox_csts*ws*ws)**(1./3.)

                if('dustvel' in key):
                    xplot=u; yplot=v; 
                    #yplot=rhod; xplot=rhog; 

                    #yplot=np.sqrt(b*b+u*u)
                    #xplot=np.sqrt(0.01*a_tilde/alpha_tilde);

                    #yplot=rhod; xplot=rhog; 
                    #xplot=rhog/mu*np.sqrt(omega_ts_kmax*omega_ts_kbox)
                    #xplot=rhog/mu*omega_ts_kbox/2.


                if('dustden' in key):
                    yplot=rhod; xplot=rhog; 

                    #xplot=np.exp(rhog)-1.; yplot=np.exp(rhod)-1.;
                    #xplot=omega_ts_kbox * Lbox_csts / 2.;
                    #xplot=np.sqrt(omega_ts_kmax*omega_ts_kbox)/0.01*xplot # sort-of-weighted average of (omega[k]/k), which is what matters
                    #yplot=np.exp(rhod)-1.; xplot=rhog; 
                    #xplot=(b/np.sqrt(3.)*np.sqrt(4.*np.pi))*(2.*np.pi/Lbox_csts)**(1./3.)*(0.01)**(-2./3.)

                    #yplot=u
                    #xplot=np.sqrt(0.01*a_tilde/alpha_tilde);



                pylab.plot(xplot,yplot,markersize=markersize,marker=marker,markerfacecolor='None',color=color,linestyle='')
                y=np.append(y,yplot); x=np.append(x,xplot);

        pylab.xscale('log'); pylab.yscale('log'); xlab=''; ylab='';
        ax=np.array([np.min(x)/2.,np.max(x)*2.,np.min(y)/2.,np.max(y)*2.])
        x0=10.**np.linspace(np.log10(ax[0]),np.log10(ax[1]),100.); y0=1.*x0;
        if('gasvel' in key):
            # some residuals here - med boxes low, small high?
            # also j very low, much weaker than this would predict; somewhat better for B-tension argument, but still low there
            # -- disperse state not very well-described by this fit
            # -- also g0,q0 very high dv versus this prediction, but unclear how significant that is
            xlab=r'Predicted Dispersion $\langle \omega \rangle\,L_{\rm box}/c_{s}$'
            ylab=r'Gas Velocity Dispersion $|\delta{\bf u}_{g}|/c_{s}$'
            xlab=r'$L_{\rm box}/(2\,c_{s}\,\langle t^{0}_{\rm grow}[L_{\rm box}] \rangle)$'
            ylab=r'$|\delta{\bf u}_{g}|/c_{s}$'
            #y0 /= 2.0;
        if('egy' in key):
            xlab=r'Kinetic Energy $\langle \rho_{g}\,|\delta{\bf u}_{g}|^{2} \rangle/2$'
            ylab=r'Magnetic Energy $\langle |\delta{\bf B}|^{2} \rangle/8\pi$'
            xlab=r'$\langle \rho_{g}\,|\delta{\bf u}_{g}|^{2} \rangle/2$'
            ylab=r'$\langle |\delta{\bf B}|^{2} \rangle/8\pi$'
        if('bmag' in key):
            xlab=r'Predicted $\rho_{d}\,|{\bf a}|\,L_{\rm box}/|{\bf B}|\,(4\pi\,P_{0})^{1/2}$'
            ylab=r'Magnetic Field $|\delta{\bf B}|/(4\pi\,P_{0})^{1/2}$'
            xlab=r'$\rho_{d}\,|{\bf a}|\,L_{\rm box}/|{\bf B}|\,(4\pi\,P_{0})^{1/2}$'
            ylab=r'$|\delta{\bf B}|/(4\pi\,P_{0})^{1/2}$'
        if('dustvel' in key):
            xlab=r'Gas Velocity Dispersion $|\delta{\bf u}_{g}|/c_{s}$'
            ylab=r'Dust Velocity Dispersion $|\delta{\bf v}_{d}|/c_{s}$'
            xlab=r'$|\delta{\bf u}_{g}|/c_{s}$'
            ylab=r'$|\delta{\bf v}_{d}|/c_{s}$'
        if('dustden' in key):
            xlab=r'Gas Density Disperion $\delta\ln(\rho_{g}/\rho_{g}^{0})$'
            ylab=r'Dust Density Disperion $\delta\ln(\rho_{d}/\rho_{d}^{0})$'
            xlab=r'$\delta\ln(\rho_{g}/\rho_{g}^{0})$'
            ylab=r'$\delta\ln(\rho_{d}/\rho_{d}^{0})$'
        if('gasden' in key):
            xlab=r'Gas Velocity Dispersion $|\delta{\bf u}_{g}|/c_{s}$'
            ylab=r'Gas Density Disperion $\delta\ln(\rho_{g}/\rho_{g}^{0})$'
            xlab=r'$|\delta{\bf u}_{g}|/c_{s}$'
            xlab=r'$|\ln\{1 + (|\delta{\bf u}_{g}|/3\,c_{s})^{2} \}|^{1/2}$'
            ylab=r'$\delta\ln(\rho_{g}/\rho_{g}^{0})$'
            #y0=np.sqrt(np.log(1.+(0.333*x0)**2.)) # 0.2-1 brackets range

        pylab.axis(ax);
        pylab.plot(x0,y0,linestyle=':',color='black',linewidth=2)
        pylab.xlabel(xlab,fontsize=f00,labelpad=2.); pylab.ylabel(ylab,fontsize=f00,labelpad=2.); 
        frame1=plot.gca() #get coordinate axes, to use for tick manipulation below
        frame1.axes.tick_params(direction='in',which='both',axis='both',
            bottom=True,top=True,left=True,right=True)
        frame1.axes.minorticks_on(); 
        frame1.axes.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10.,numticks=55))
        frame1.axes.tick_params(axis='y',pad=1)

    pylab.subplots_adjust(left=0.1,bottom=0.1,right=1-0.01,top=1-0.01,hspace=0.2,wspace=0.25);
    pylab.savefig('dust_corr.pdf',transparent=True,bbox_inches='tight',pad_inches=0)



def crunchalldispersions():
    s0='/panfs/ds08/hopkins/phopkins/data2/dust_instability_boxes/mhd_rhi_dust_boxes/'
    sdirs=['AGB/x0_N128','AGB/u0_N128','AGB/v0_N128','AGB/w0_N128_fixed','CGM_QSO_Lkpc/j0_N128_boris','CGM_QSO_Lkpc/j0_N128_gamma53_pelecscalecharge_boris',
    'HII_inner_L0pt1pc/k0_N128','HII_inner_L0pt1pc/k0_N128_pelecscalecharge','HII_inner_L0pt1pc/k0_N128_pelecscalecharge_2',
    'HII_inner_Lau/c0_N128','HII_inner_Lsubau/t0_N128_fixed','HII_outer_L200au/s0_N128','HII_outer_Lpc/e0_N128',
    'jono_default_params/b0_N128','jono_default_params/b0_N128_beta1traditionalbeta','jono_default_params/b0_N128_gamma53',
    'WIM_Lpc/f0_N128_boris','WIM_Lpc/f0_N128_5x','WIM_L500au/q0_N128_boris',
    'HII_inner_L0pt1pc/k0_N128_mu','HII_inner_L0pt1pc/k0_N128_mu2',
    'jono_default_params/b0_N128_mu','jono_default_params/b0_N128_mu2']

    sdirs=['active_runs/g0_N128_5x','active_runs/m0_N128','active_runs/g0_N128','active_runs/q0_N128_5x']
    sdirs=['HII_inner_L0pt1pc/k0_N128_betacorr']
    sdirs=['CGM_QSO_Lkpc/j0_N128_betacorr']

    for sdir in sdirs:
        crunchdispersions(sdir=sdir,sdir_master=s0)


def crunchdispersions(sdir,n_touse=5,sdir_master='/panfs/ds08/hopkins/phopkins/data2/dust_instability_boxes/mhd_rhi_dust_boxes/',sdir_sub='/output/'):
    sdir=sdir_master+sdir+sdir_sub
    snums=np.array([],dtype='int')
    for snum in range(800):
        fname,fname_base,fname_ext = check_if_filename_exists(sdir,snum)
        if(fname != 'NULL'): snums=np.append(snums,snum)
    if(snums.size < 10*n_touse): n_touse=snums.size/10
    if(n_touse < 1): n_touse=1;
    snums=snums[snums.size-n_touse::]
    dv_g_vwt=np.zeros(3); dv_g_mwt=np.zeros(3); dv_d_vwt=np.zeros(3); dv_d_mwt=np.zeros(3); dB_g_vwt=np.zeros(3); dB_g_mwt=np.zeros(3);
    drho_g_vwt=0.; drho_g_mwt=0.; drho_d_vwt=0.; drho_d_mwt=0.; drho_dg_vwt=0.; drho_dg_mwt=0.;
    print(snums)
    for snum in snums:
        print('scanning snapshot ',snum)
        F=load_snap(sdir,snum);
        for ptype in ['PartType0','PartType3']:
            P=F[ptype]; V=P['Velocities']; vx=V[:,0]; vy=V[:,1]; vz=V[:,2]; m=np.array(P['Masses']);
            if(ptype=='PartType0'):
                rho=np.array(P['Density']); B=P['MagneticField']; bx=B[:,0]; by=B[:,1]; bz=B[:,2];
            else:
                rho,h_d,rho_gn,h_gn=load_snap_dustrho(sdir,snum)
            rho += 1.e-10; m[(m <= 0)] += m[(m>0.)].min(); rho[np.isfinite(rho)==False]=rho[np.isfinite(rho)==True].min()
            vol=m/rho; wt_v=vol/np.sum(vol); wt_m=m/np.sum(m); lnrho=np.log(rho);
            if(ptype=='PartType0'):
                dv_g_vwt[0]+=std_wt(vx,wt_v); dv_g_vwt[1]+=std_wt(vy,wt_v); dv_g_vwt[2]+=std_wt(vz,wt_v);
                dv_g_mwt[0]+=std_wt(vx,wt_m); dv_g_mwt[1]+=std_wt(vy,wt_m); dv_g_mwt[2]+=std_wt(vz,wt_m);
                dB_g_vwt[0]+=std_wt(bx,wt_v); dB_g_vwt[1]+=std_wt(by,wt_v); dB_g_vwt[2]+=std_wt(bz,wt_v);
                dB_g_mwt[0]+=std_wt(bx,wt_m); dB_g_mwt[1]+=std_wt(by,wt_m); dB_g_mwt[2]+=std_wt(bz,wt_m);
                drho_g_vwt+=std_wt(lnrho,wt_v); drho_g_mwt+=std_wt(lnrho,wt_m);
            else:
                dv_d_vwt[0]+=std_wt(vx,wt_v); dv_d_vwt[1]+=std_wt(vy,wt_v); dv_d_vwt[2]+=std_wt(vz,wt_v);
                dv_d_mwt[0]+=std_wt(vx,wt_m); dv_d_mwt[1]+=std_wt(vy,wt_m); dv_d_mwt[2]+=std_wt(vz,wt_m);
                drho_d_vwt+=std_wt(lnrho,wt_v); drho_d_mwt+=std_wt(lnrho,wt_m);
                lnrho_dg=np.log(rho/rho_gn);
                drho_dg_vwt+=std_wt(lnrho_dg,wt_v); drho_dg_mwt+=std_wt(lnrho_dg,wt_m);
        F.close()
    print(dv_d_mwt)
    print(dv_d_mwt.min())
    print(dv_d_mwt.max())
    ninv=1./(1.*snums.size);
    dv_g_vwt=np.sqrt(dv_g_vwt*ninv); dv_g_mwt=np.sqrt(dv_g_mwt*ninv);
    dv_d_vwt=np.sqrt(dv_d_vwt*ninv); dv_d_mwt=np.sqrt(dv_d_mwt*ninv);
    dB_g_vwt=np.sqrt(dB_g_vwt*ninv); dB_g_mwt=np.sqrt(dB_g_mwt*ninv);
    drho_g_vwt=np.sqrt(drho_g_vwt*ninv); drho_g_mwt=np.sqrt(drho_g_mwt*ninv);
    drho_d_vwt=np.sqrt(drho_d_vwt*ninv); drho_d_mwt=np.sqrt(drho_d_mwt*ninv);
    drho_dg_vwt=np.sqrt(drho_dg_vwt*ninv); drho_dg_mwt=np.sqrt(drho_dg_mwt*ninv);
    print('SDIR == ',sdir)
    print('delta_v_gas[V-wt] = ',dv_g_vwt)
    print('delta_v_gas[M-wt] = ',dv_g_mwt)
    print('delta_B_gas[V-wt] = ',dB_g_vwt)
    print('delta_B_gas[M-wt] = ',dB_g_mwt)
    print('delta_lnrho_gas[V-wt] = ',drho_g_vwt)
    print('delta_lnrho_gas[M-wt] = ',drho_g_mwt)
    print('delta_v_dust[V-wt] = ',dv_d_vwt)
    print('delta_v_dust[M-wt] = ',dv_d_mwt)
    print('delta_lnrho_dust[V-wt] = ',drho_d_vwt)
    print('delta_lnrho_dust[M-wt] = ',drho_d_mwt)
    print('delta_lnrho_dust/gas[V-wt] = ',drho_dg_vwt)
    print('delta_lnrho_dust/gas[M-wt] = ',drho_dg_mwt)

            
def std_wt(x,wt):
    ok=np.where( (np.isfinite(x)==True) & (np.isnan(x)==False) & (np.isfinite(wt)==True) & (np.isnan(wt)==False) & (np.abs(x) < 1.e10))
    w=1.*wt[ok]/np.sum(1.*wt[ok]); y=1.*x[ok];
    ok=np.where( (np.isfinite(y)==True) & (np.isnan(y)==False) & (np.isfinite(w)==True) & (np.isnan(w)==False)  & (np.abs(y) < 1.e10))
    w=1.*w[ok]; y=1.*y[ok]; f=1.*y*w;
    q=np.sum(f); 
    return np.sum(y*f)-q*q          


def return_tsi(sdir):
    tsi=0.;
    if('b0_' in sdir): tsi=0.34
    if('w0_' in sdir): tsi=0.01
    if('v0_' in sdir): tsi=3.1
    if('u0_' in sdir): tsi=920.
    if('x0_' in sdir): tsi=3.e5
    if('t0_' in sdir): tsi=0.0088
    if('c0_' in sdir): tsi=2.6
    if('k0_' in sdir): tsi=1000.
    if('l0_' in sdir): tsi=0.0032
    if('s0_' in sdir): tsi=1.0
    if('e0_' in sdir): tsi=330.
    if('g0_' in sdir): tsi=4.8e-4
    if('q0_' in sdir): tsi=0.21
    if('f0_' in sdir): tsi=100.
    if('m0_' in sdir): tsi=3.3e-4
    if('r0_' in sdir): tsi=0.067
    if('h0_' in sdir): tsi=20.
    if('j0_' in sdir): tsi=0.29
    return tsi


def crunchfigs(gas_val_toplot='b',
        sdir='/panfs/ds08/hopkins/phopkins/data2/dust_instability_boxes/mhd_rhi_dust_boxes/HII_outer_Lpc/e0_N128/output',
        vmin=0,vmax=0,alpha=0.1,markersize=5.,
        snums=range(800),resolution='high',cmap='terrain'):
    for snum in snums:
        fname,fname_base,fname_ext = check_if_filename_exists(sdir,snum)
        if(fname != 'NULL'):
            plotpts_3d_2d(snum=snum,sdir=sdir,forsavedfigure=True,vmin=vmin,vmax=vmax,
                resolution=resolution,cmap=cmap,gas_val_toplot=gas_val_toplot,
                ts0_inv=return_tsi(sdir),alpha=alpha,markersize=markersize)
        


def plotpts_3d_2d(snum=0,sdir='./output',ptype='PartType3',width=0.05,cut_dust=1.,
        vmin=0,vmax=0,forsavedfigure=False,gas_val_toplot='b',zmed_set=0.,alpha=0.1,markersize=5.,
        max_ngas_tokeep=1.e6,resolution='low',cmap='terrain',fontsize=26.,
        ts0_inv=0.,center_acc_vec=np.zeros(0),cen_list=np.zeros(0),time_list=np.zeros(0)):
    pylab.close('all'); fig=plot.figure(1,figsize=(12.,12.)); ax=fig.add_subplot(111,projection='3d');
    ax.set_xlim3d(0,1); ax.set_ylim3d(0,1); ax.set_zlim3d(0,1);

    P_File=load_snap(sdir,snum);
    time_snap = P_File['Header'].attrs['Time']
    if(ts0_inv>0.): time_snap *= ts0_inv
    print('Time == ',time_snap)
    tstr = "%.0f"%time_snap
    if(time_snap < 10.): tstr = "%.1f"%time_snap
    if(time_snap < 0.1): tstr = "%.2f"%time_snap
    if(time_snap < 0.01): tstr = "%.3f"%time_snap
    if(time_snap < 0.001): tstr = "%.4f"%time_snap
    if(time_snap < 0.0001): tstr = "%.5f"%time_snap
    if(ts0_inv>0.): tstr = tstr + r'\,$t_{s}^{0}$'
    ax.text2D(0.08,0.96,r'${\rm Time}=\,$'+tstr,color='black',
        fontsize=fontsize*1.15,ha='left',va='top',transform=ax.transAxes)
    P = P_File[ptype]
    print('Dust-to-Gas Mass Ratio = ',np.sum(P_File['PartType3']['Masses'][:])/np.sum(P_File['PartType0']['Masses'][:]))
    Pc=np.array(P['Coordinates']);
    center=np.zeros(0); center_acc_vec=np.array(center_acc_vec)
    if(center_acc_vec.size > 0):
        time_code = P_File['Header'].attrs['Time']
        center = 0.5 + 0.5 * center_acc_vec * time_code*time_code
        center = center - np.floor(center)
    if(time_list.size > 0):
        time_code = P_File['Header'].attrs['Time']
        center=np.zeros(3); 
        center[0]=np.interp(time_code,time_list,cen_list[:,0])
        center[1]=np.interp(time_code,time_list,cen_list[:,1])
        center[2]=np.interp(time_code,time_list,cen_list[:,2])
        center = center - np.floor(center)
    if(center.size > 0): Pc-=center;
    Pc[(Pc>1.)] -= 1.;
    Pc[(Pc<0.)] += 1.;
    vx = np.array(P['Velocities'][:,0])
    vy = np.array(P['Velocities'][:,1])
    vz = np.array(P['Velocities'][:,2])
    print('Var v_x_dust == ',np.min(vx),np.median(vx),np.max(vx),np.std(vx),' (min/median/max/std)')
    print('Var v_y_dust == ',np.min(vy),np.median(vy),np.max(vy),np.std(vy),' (min/median/max/std)')
    print('Var v_z_dust == ',np.min(vz),np.median(vz),np.max(vz),np.std(vz),' (min/median/max/std)')
    vmin_prev=vmin; vmax_prev=vmax;
    for subplot,xz,yz,xyzpanel,zpostoplot in zip([1,2,3],[0,1,0],[0,0,1],['z','y','x'],[1.,0.,1.]): 
        xplotc=0; yplotc=1; depth_c=2; #plot xy, z=depth
        if(xz==1): 
            xplotc=0; yplotc=2; depth_c=1; #plot xz, y=depth
        if(yz==1): 
            xplotc=1; yplotc=2; depth_c=0; #plot yz, x=depth
        
        quiet=True; 
        if(subplot==1): quiet=False # give numbers once
        z=Pc[:,depth_c]; zmx=np.max(z)-np.min(z); z0=np.median(z); 
        if(zmed_set > -1.e9): z0=zmed_set 
        zz=np.abs(z-z0); zz[(zz>0.5*zmx)] = zmx-zz[(zz>0.5*zmx)]
        ok_master = (np.random.rand(z.size) < (1.*max_ngas_tokeep)/(1.+1.*z.size))
        ok=np.where((zz < width) & ok_master)
        x=Pc[:,xplotc][ok]; y=Pc[:,yplotc][ok]; 
        ok_r = (np.random.rand(x.size) < cut_dust)
        vmin_0=vmin; vmax_0=vmax; 
        if(vmin_0==0): 
            if(subplot>1): vmin=vmin_prev
        if(vmax_0==0):
            if(subplot>1): vmax=vmax_prev
        vmin_prev, vmax_prev, im_save = gas_rho_image_multidim(ax, snum=snum,sdir=sdir,xmax=1.,
            xz=xz,yz=yz,gas_val_toplot=gas_val_toplot,zmed_set=z0,fontsize=fontsize,
            vmin=vmin,vmax=vmax,quiet=quiet,center=center,
            zdir=xyzpanel,zs=zpostoplot,resolution=resolution,cmap=cmap)
        ax.plot(x[ok_r],y[ok_r],marker='.',markersize=markersize,alpha=alpha,
            linestyle='',color='black',markeredgecolor='None',zdir=xyzpanel,zs=zpostoplot)
        ax.xaxis.set_ticklabels([]); ax.yaxis.set_ticklabels([]); #no tick names
        ax.get_xaxis().set_ticks([]); ax.get_yaxis().set_ticks([]); ax.set_zticks([]); #no ticks
    ax.dist=9
    ext='000'+str(snum);
    if (snum>=10): ext='00'+str(snum)
    if (snum>=100): ext='0'+str(snum)
    if (snum>=1000): ext=str(snum)
    fname=sdir+'/'+'im_3d2dProj_'+gas_val_toplot+'_'+resolution+'_'+cmap+'_'+ext
    if(forsavedfigure==True): 
        fig.savefig(fname+'.png',dpi=150,bbox_inches='tight',pad_inches=0)
    P_File.close()




def plotpts_w_gas_loop(snum_min=0,snum_max=3300,sdir='./output',width=0.05,
        cut_dust=0.75,vmin=-1.,vmax=1.,fontsize=20.):
    for snum in np.arange(snum_min,snum_max+1):
        plotpts_w_gas(snum=snum,sdir=sdir,cut_dust=cut_dust,width=width,vmin=vmin,vmax=vmax)


def crunchfigs_movie(snums=range(700),resolution='high',cmap='terrain',cmap_dust='inferno',
        #sdir='HII_inner_Lau/c0_N128',gas_val_toplot='b',vmin=0.05,vmax=0.75,vmin_dust=-1.5,vmax_dust=6.2,center_acc_vec=[10.6066,0.,10.6066],alpha=0.1,markersize=5.,max_ngas_tokeep=1e6,
        #sdir='jono_default_params/b0_N128',gas_val_toplot='v',vmin=0.001,vmax=0.22,vmin_dust=-1.5,vmax_dust=5.,center_acc_vec=[4.99315,0.,0.26168],alpha=0.1,markersize=5.,max_ngas_tokeep=1e6,
        #sdir='HII_outer_Lpc/e0_N128',gas_val_toplot='b',vmin=0.1,vmax=1.5,vmin_dust=-1.5,vmax_dust=10.,center_acc_vec=[63.63961,0.,63.63961],alpha=0.1,markersize=5.,max_ngas_tokeep=1e6,
        sdir='HII_inner_L0pt1pc/k0_MG_N128_Ac',gas_val_toplot='rho',vmin=-5.,vmax=5.,vmin_dust=-1.5,vmax_dust=10.,center_acc_vec=[63.63961,0.,63.63961],alpha=0.15,markersize=1.5,max_ngas_tokeep=1e8,
        sdir_master='/panfs/ds08/hopkins/phopkins/data2/dust_instability_boxes/mhd_rhi_dust_boxes/'):
        
    sdir_full = sdir_master + sdir + '/output/'
    
    ## first we need to do a loop to get the centering right:
    time_list=np.zeros(0); cen_list=np.zeros((0,3)); cen=np.array([0.5,0.5,0.5]); vel_prev=np.zeros(3); vel_new=np.zeros(3); t_prev=0.; t_new=0.;
    fname_cenlist = sdir_full + 'movie_cen_save.hdf5'
    if not os.path.exists(fname_cenlist):
        print("No centroid-list exists: generating...")
        for snum in range(700):
            fname,fname_base,fname_ext = check_if_filename_exists(sdir_full,snum)
            if(fname != 'NULL'):
                ptype='PartType0'
                F=load_snap(sdir_full,snum); t_prev=t_new; t_new=F['Header'].attrs['Time'];
                m=np.array(F[ptype]['Masses']); p=np.array(F[ptype]['Density']); v=np.array(F[ptype]['Velocities']);
                wt=m/p; wt[np.isfinite(wt)==False]=0.; wt/=np.sum(wt); vel_prev=vel_new; vel_new=np.array([np.sum(wt*v[:,0]),np.sum(wt*v[:,1]),np.sum(wt*v[:,2])])
                cen+=0.5*(vel_new+vel_prev)*(t_new-t_prev);
                time_list=np.append(time_list,t_new)
                cen_list=np.concatenate((cen_list,[cen]),axis=0)
        file = h5py.File(fname_cenlist,'w') 
        file.create_dataset("times",data=time_list)
        file.create_dataset("centroids",data=cen_list)
        file.close()
    else:
        print("Pre-compiled centroid-list exists: using...")
        file = h5py.File(fname_cenlist,'r') 
        time_list = 1.*np.array(file['times'])
        cen_list = 1.*np.array(file['centroids'])
        file.close()
    
    for snum in snums:
        fname,fname_base,fname_ext = check_if_filename_exists(sdir_full,snum)
        if(fname != 'NULL'):
            print('crunching snapshot ',snum)
            if(2==2):
                plotpts_3d_2d(snum=snum,sdir=sdir_full,forsavedfigure=True,vmin=vmin,vmax=vmax,max_ngas_tokeep=max_ngas_tokeep,
                    time_list=time_list,cen_list=cen_list,center_acc_vec=center_acc_vec,alpha=alpha,markersize=markersize,
                    resolution=resolution,cmap=cmap,gas_val_toplot=gas_val_toplot,ts0_inv=return_tsi(sdir))
            if(2==0):
                plotpts_3dscatter(snum=snum,sdir=sdir_full,vmin=vmin_dust,vmax=vmax_dust,cmap=cmap_dust,
                    ts0_inv=return_tsi(sdir),time_list=time_list,cen_list=cen_list,center_acc_vec=center_acc_vec,
                    forsavedfigure=True)
        



def plotpts_3dscatter(snum=0,sdir='./output/jonodefault_b0_128',ptype='PartType3',
        gas_val_toplot='rho',number_to_plot=1.e8,alpha=0.05,
        vmin=0,vmax=0,fontsize=26.,forsavedfigure=True,time_list=np.zeros(0),cen_list=np.zeros(0),
        cmap='inferno',rasterized=True,ts0_inv=0.,center_acc_vec=np.zeros(0)):
        #-3,2 vs -1.87,5.1

    pylab.close('all'); fig=plot.figure(1,figsize=(12.,12.)); ax=fig.add_subplot(111,projection='3d');
    P_File=load_snap(sdir,snum);
    time_snap = P_File['Header'].attrs['Time']
    if(ts0_inv>0.): time_snap *= ts0_inv
    print('Time == ',time_snap)
    tstr = "%.0f"%time_snap
    if(time_snap < 10.): tstr = "%.1f"%time_snap
    if(time_snap < 0.1): tstr = "%.2f"%time_snap
    if(time_snap < 0.01): tstr = "%.3f"%time_snap
    if(time_snap < 0.001): tstr = "%.4f"%time_snap
    if(time_snap < 0.0001): tstr = "%.5f"%time_snap
    if(ts0_inv>0.): tstr = tstr + r'\,$t_{s}^{0}$'
    ax.text2D(0.08,0.96,r'${\rm Time}=\,$'+tstr,color='black',
        fontsize=fontsize*1.15,ha='left',va='top',transform=ax.transAxes)
    P = P_File[ptype]; Pc=np.array(P['Coordinates']); Pv=np.array(P['Velocities']);
    center=np.zeros(0); center_acc_vec=np.array(center_acc_vec); time_list=np.array(time_list); cen_list=np.array(cen_list);
    if(center_acc_vec.size > 0):
        time_code = P_File['Header'].attrs['Time']
        center = 0.5 + 0.5 * center_acc_vec * time_code*time_code
        center = center - np.floor(center)
    if(time_list.size > 0):
        time_code = P_File['Header'].attrs['Time']
        center=np.zeros(3); 
        center[0]=np.interp(time_code,time_list,cen_list[:,0])
        center[1]=np.interp(time_code,time_list,cen_list[:,1])
        center[2]=np.interp(time_code,time_list,cen_list[:,2])
        center = center - np.floor(center)
    if(center.size > 0): Pc-=center;
    Pc[(Pc>1.)] -= 1.;
    Pc[(Pc<0.)] += 1.;
    x=Pc[:,0]; y=Pc[:,1]; z=Pc[:,2]; vx=Pv[:,0]; vy=Pv[:,1]; vz=Pv[:,2];
    dvx=vx-np.median(vx); dvy=vy-np.median(vy); dvz=vz-np.median(vz); v=np.sqrt(dvx*dvx+dvy*dvy+dvz*dvz); 
    ok=np.where(np.random.rand(x.size) < number_to_plot/(1.+x.size))
    colorv=v
    sv=5.
    if(gas_val_toplot=='rho'): 
        if(ptype=='PartType0'): 
            colorv=np.array(P['Density'])
        else:
            rho_d, h_d, rho_gn, h_gn = load_snap_dustrho(sdir,snum)
            colorv = np.log(rho_d / 0.01)
    print('min_max range of colorv = ',np.min(colorv),np.median(colorv),np.max(colorv))
    s = np.argsort(colorv)
    if((vmin==0)|(vmax==0)):
        vmed=colorv[s][np.round(0.5*s.size).astype('int')]
        if(vmin==0): 
            vmin=colorv[s][np.round(0.01*s.size).astype('int')]
            vmm=vmed-0.001;
            if(vmin>vmm): vmin=vmm
        if(vmax==0): 
            vmax=colorv[s][np.round(0.999*s.size).astype('int')]
            vmm=vmed+0.001;
            if(vmax<vmm): vmax=vmm
        
    colorv[(colorv < vmin)]=vmin; colorv[(colorv > vmax)]=vmax;
    ax.set_xlim3d(0,1); ax.set_ylim3d(0,1); ax.set_zlim3d(0,1);
    im_for_cbar = ax.scatter(x[ok],y[ok],z[ok],
        marker='.',edgecolor='face',s=10.,
        rasterized=rasterized,c=colorv[ok],cmap=cmap,
        alpha=alpha,vmin=vmin,vmax=vmax)
        
    ax.xaxis.set_ticklabels([]); ax.yaxis.set_ticklabels([]); #no tick names
    ax.get_xaxis().set_ticks([]); ax.get_yaxis().set_ticks([]); ax.set_zticks([]); #no ticks
    labelcbar = r'$\ln{[\rho_{\rm d}/\rho_{\rm d}^{0}]}$'

    norm=matplotlib.colors.Normalize(vmin=vmin,vmax=vmax)
    m = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    m.set_array(colorv)
    cb1=pylab.colorbar(m,shrink=0.7,aspect=30,fraction=0.05,pad=-0.04)
    cb1.set_label(labelcbar,fontsize=fontsize*1.15,labelpad=10.)
    cb1.ax.tick_params(labelsize=fontsize)
    cb1.ax.yaxis.set_tick_params(direction='in')
    ax.xaxis.set_pane_color((0.0, 0.0, 0.0, 0.5))
    ax.yaxis.set_pane_color((0.0, 0.0, 0.0, 0.5))
    ax.zaxis.set_pane_color((0.0, 0.0, 0.0, 0.5))

    ax.dist=9
    ext='000'+str(snum);
    if (snum>=10): ext='00'+str(snum)
    if (snum>=100): ext='0'+str(snum)
    if (snum>=1000): ext=str(snum)
    fname=sdir+'/'+'im_3dDust_'+cmap+'_'+ext
    if(forsavedfigure==True): 
        pylab.savefig(fname+'.png',dpi=150,bbox_inches='tight',pad_inches=0)
    P_File.close()


def plotpts_3d_render(snum=0,sdir='./output/',ptype='PartType3',
        gas_val_toplot='rho',number_to_plot=1.e8,alpha=0.1,
        vmin=-3.,vmax=2.,fontsize=26.,forsavedfigure=False,
        cmap='terrain',gridpt_num=256j,
        opacity=0.2,ncontours=8):

    pylab.close('all');
    #fig=plot.figure(1,figsize=(12.,12.)); 
    #ax=fig.add_subplot(111,projection='3d');
    P_File=load_snap(sdir,snum);
    time_snap = P_File['Header'].attrs['Time']
    print('Time == ',time_snap)
    tstr = "%.2f"%time_snap
    if(time_snap < 0.01): tstr = "%.3f"%time_snap
    if(time_snap < 0.001): tstr = "%.4f"%time_snap
    if(time_snap < 0.0001): tstr = "%.5f"%time_snap
    #ax.text2D(0.08,0.96,r'${\rm Time}=\,$'+tstr,color='black',fontsize=fontsize,ha='left',va='top',transform=ax.transAxes)
    P = P_File[ptype]; Pc=np.array(P['Coordinates']); Pv=np.array(P['Velocities']);
    x=Pc[:,0]; y=Pc[:,1]; z=Pc[:,2]; vx=Pv[:,0]; vy=Pv[:,1]; vz=Pv[:,2];
    dvx=vx-np.median(vx); dvy=vy-np.median(vy); dvz=vz-np.median(vz); v=np.sqrt(dvx*dvx+dvy*dvy+dvz*dvz); 
    ok=np.where(np.random.rand(x.size) < number_to_plot/(1.+x.size))
    colorv=v
    sv=5.
    if(gas_val_toplot=='rho'): 
        if(ptype=='PartType0'): 
            colorv=np.array(P['Density'])
        else:
            rho_d, h_d, rho_gn, h_gn = load_snap_dustrho(sdir,snum)
            colorv = np.log(rho_d / 0.01)
    print('min_max range of colorv = ',np.min(colorv),np.median(colorv),np.max(colorv))
    xg,yg,zg = np.mgrid[0:1:gridpt_num, 0:1:gridpt_num, 0:1:gridpt_num] # create grid
    print('entering interp')
    F = interpolate.griddata(Pc,colorv,(xg,yg,zg),fill_value=np.min(colorv)-2.5,method='nearest') # interp to grid
    print('done with interp')
    colorv=F; #colorv[(colorv < vmin)]=vmin; colorv[(colorv > vmax)]=vmax; # clip the data
    # accent, 

    import mayavi.mlab as mlab
    from tvtk.util import ctf
    #mlab.contour3d(colorv,contours=8,opacity=.2,colormap=cmap,vmin=vmin,vmax=vmax)
    if(2==0):
        mlab.contour3d(colorv,contours=ncontours,opacity=opacity,
            colormap=cmap,vmin=vmin,vmax=vmax)
    else:
        
        scalar_data = mlab.pipeline.scalar_field(colorv)
        
        if(2==2):
            volume = mlab.pipeline.volume(scalar_data,vmin=vmin,vmax=vmax)
            
            #values = np.linspace(0., 1., 256) # dummy field for map conversion
            #c = ctf.save_ctfs(volume._volume_property) # save the existing colormap
            #c['rgb']=plot.cm.get_cmap(cmap)(values.copy()) # change it with the colors of the new colormap
            #ctf.load_ctfs(c, volume._volume_property) # load the color transfer function to the volume
            #volume.update_ctf = True # signal for update
        else:
            cval_v = np.linspace(vmin,vmax,ncontours)
            oval_v = np.linspace(0.01,1.,ncontours)
            oval_v = np.linspace(0.05,0.5,ncontours)
            oval_v *= oval_v
            for cv,ov in zip(cval_v,oval_v):        
                mlab.pipeline.iso_surface(scalar_data,vmin=vmin,vmax=vmax,
                    colormap=cmap,contours=[cv],opacity=ov)
        #    contours=[s.min()+0.1*s.ptp(), ], opacity=0.3)
        #mlab.pipeline.iso_surface(src, contours=[s.max()-0.1*s.ptp(), ],)
        
        # below block is just to change the colormaps
        
    #ax.set_xlim3d(0,1); ax.set_ylim3d(0,1); ax.set_zlim3d(0,1);

    print('done-ish')


def tester():
    # Create some test data, 3D gaussian, 200 points
    dx, pts = 2, 100j

    N = 500
    R = np.random.random((N,3))*2*dx - dx
    V = np.exp(-( (R**2).sum(axis=1)) )

    # Create the grid to interpolate on
    X,Y,Z = np.mgrid[-dx:dx:pts, -dx:dx:pts, -dx:dx:pts]

    # Interpolate the data
    F = griddata(R, V, (X,Y,Z))
    
    # From here it's a snap to display our data:
    contour3d(F,contours=8,opacity=.2,colormap='spectral')


    ax.xaxis.set_ticklabels([]); ax.yaxis.set_ticklabels([]); #no tick names
    ax.get_xaxis().set_ticks([]); ax.get_yaxis().set_ticks([]); ax.set_zticks([]); #no ticks
    labelcbar = r'$\ln{[\rho_{\rm d}/\rho_{\rm d}^{0}]}$'

    norm=matplotlib.colors.Normalize(vmin=vmin,vmax=vmax)
    m = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    m.set_array(colorv)
    cb1=pylab.colorbar(m,shrink=0.7,aspect=30,fraction=0.05,pad=-0.04)
    cb1.set_label(labelcbar,fontsize=fontsize,labelpad=10.)
    cb1.ax.tick_params(labelsize=fontsize)

    ax.dist=9
    ext='000'+str(snum);
    if (snum>=10): ext='00'+str(snum)
    if (snum>=100): ext='0'+str(snum)
    if (snum>=1000): ext=str(snum)
    fname=sdir+'/'+'im_3dDust_'+cmap+'_'+ext
    if(forsavedfigure==True): pylab.savefig(fname+'.png',dpi=150,bbox_inches='tight',pad_inches=0)
    P_File.close()



def plotpts(snum=0,sdir='./output',ptype='PartType0',width=0.05):
    pylab.close('all'); 
    plot.figure(1,figsize=(16.,8.))
    P_File=load_snap(sdir,snum);
    print('Time == ',P_File['Header'].attrs['Time'])
    P = P_File[ptype]
    Pc=np.array(P['Coordinates']);
    z=Pc[:,2]; z0=np.median(z);
    zz=np.abs(z-z0)
    ok=np.where(zz < width)
    x=Pc[:,0][ok]
    y=Pc[:,1][ok]
    z=Pc[:,2][ok]
    pylab.subplot(1,2,1)
    pylab.plot(x,y,marker=',',linestyle='',rasterized=True)
    z=Pc[:,1]; z0=np.median(z);
    zz=np.abs(z-z0)
    ok=np.where(zz < width)
    x=Pc[:,0][ok]
    y=Pc[:,1][ok]
    z=Pc[:,2][ok]
    pylab.subplot(1,2,2)
    pylab.plot(x,z,marker=',',linestyle='',rasterized=True)
    vx = P['Velocities'][:,0]
    vy = P['Velocities'][:,1]
    vz = P['Velocities'][:,2]
    print('Var v_x == ',np.min(vx),np.median(vx),np.max(vx),np.std(vx))
    print('Var v_y == ',np.min(vy),np.median(vy),np.max(vy),np.std(vy))
    print('Var v_z == ',np.min(vz),np.median(vz),np.max(vz),np.std(vz))

    vy_g = P_File['PartType0']['Velocities'][:,1]
    print('Mean v_y_gas == ',np.mean(vy_g))
    print('v_drift == ',np.mean(vy)-np.mean(vy_g),np.median(vy)-np.median(vy_g))

    if(ptype=='PartType0'):
        rho = P['Density'][:]
        print('Densities == ',np.min(rho),np.median(rho),np.max(rho),np.std(rho))
    P_File.close()



def load_snap(sdir,snum,snapshot_name='snapshot',extension='.hdf5',four_char=0):
    fname,fname_base,fname_ext = check_if_filename_exists(sdir,snum,\
        snapshot_name=snapshot_name,extension=extension,four_char=four_char)
    file = h5py.File(fname,'r')
    return file


def load_snap_dustrho(sdir,snum,snapshot_name='dust_snap',extension='.hdf5',four_char=1):
    fname,fname_base,fname_ext = check_if_filename_exists(sdir,snum,\
        snapshot_name=snapshot_name,extension=extension,four_char=four_char)
    if(fname == 'NULL'):
        import grains.grain_density_from_snapshot as gd
        h_d,nngb_d,h_g,nngb_g=gd.grain_density_from_snapshot(sdir=sdir,snum=snum,
            Hmin=0.,boxsize_x=1,boxsize_y=1,boxsize_z=1)
        fname,fname_base,fname_ext = check_if_filename_exists(sdir,snum,\
            snapshot_name=snapshot_name,extension=extension,four_char=four_char)
    if(fname != 'NULL'):
        file = h5py.File(fname,'r')
        q = file['Mass_Density_List_Of_Dust_Neighbors'][:]
        h = file['Smoothing_Length_List_Of_Dust_Neighbors'][:]
        qg = file['Mass_Density_List_Of_Gas_Neighbors'][:]
        hg = file['Smoothing_Length_List_Of_Gas_Neighbors'][:]
        file.close()
        return q, h, qg, hg



def check_if_filename_exists(sdir,snum,snapshot_name='snapshot',extension='.hdf5',four_char=0):
    for extension_touse in [extension,'.bin','']:
        fname=sdir+'/'+snapshot_name+'_'
        ext='00'+str(snum);
        if (snum>=10): ext='0'+str(snum)
        if (snum>=100): ext=str(snum)
        if (four_char==1): ext='0'+ext
        if (snum>=1000): ext=str(snum)
        fname+=ext
        fname_base=fname

        s0=sdir.split("/"); snapdir_specific=s0[len(s0)-1];
        if(len(snapdir_specific)<=1): snapdir_specific=s0[len(s0)-2];

        ## try several common notations for the directory/filename structure
        fname=fname_base+extension_touse;
        if not os.path.exists(fname): 
            ## is it a multi-part file?
            fname=fname_base+'.0'+extension_touse;
        if not os.path.exists(fname): 
            ## is the filename 'snap' instead of 'snapshot'?
            fname_base=sdir+'/snap_'+ext; 
            fname=fname_base+extension_touse;
        if not os.path.exists(fname): 
            ## is the filename 'snap' instead of 'snapshot', AND its a multi-part file?
            fname=fname_base+'.0'+extension_touse;
        if not os.path.exists(fname): 
            ## is the filename 'snap(snapdir)' instead of 'snapshot'?
            fname_base=sdir+'/snap_'+snapdir_specific+'_'+ext; 
            fname=fname_base+extension_touse;
        if not os.path.exists(fname): 
            ## is the filename 'snap' instead of 'snapshot', AND its a multi-part file?
            fname=fname_base+'.0'+extension_touse;
        if not os.path.exists(fname): 
            ## is it in a snapshot sub-directory? (we assume this means multi-part files)
            fname_base=sdir+'/snapdir_'+ext+'/'+snapshot_name+'_'+ext; 
            fname=fname_base+'.0'+extension_touse;
        if not os.path.exists(fname): 
            ## is it in a snapshot sub-directory AND named 'snap' instead of 'snapshot'?
            fname_base=sdir+'/snapdir_'+ext+'/'+'snap_'+ext; 
            fname=fname_base+'.0'+extension_touse;
        if not os.path.exists(fname): 
            ## wow, still couldn't find it... ok, i'm going to give up!
            fname_found = 'NULL'
            fname_base_found = 'NULL'
            fname_ext = 'NULL'
            continue;
        fname_found = fname;
        fname_base_found = fname_base;
        fname_ext = extension_touse
        break; # filename does exist! 
    return fname_found, fname_base_found, fname_ext;




def gas_rho_image(snum=3,sdir='./output',vmin=0.,vmax=0.,ptype='PartType0',
        cmap='terrain',xmax=1.,xz=0,yz=0,gas_val_toplot='rho',rhocut=-2.5,save='dummy',
        zmed_set=-1.e10,quiet=False,zdir='z',zs=0.):
    P_File=load_snap(sdir,snum);
    P = P_File[ptype]
    Pc= np.array(P['Coordinates'])
    dmax = np.max(Pc)
    xset=0; yset=1; zset=2;
    if(xz==1): 
        xset=0; yset=2; zset=1;
    if(yz==1): 
        xset=1; yset=2; zset=0;
    xx = Pc[:,xset]; yy = Pc[:,yset]; zz = Pc[:,zset]
    vx = P['Velocities'][:,0]; vy = P['Velocities'][:,1]; vz = P['Velocities'][:,2]
    if(quiet==False):
        print('Var v_x_gas == ',np.min(vx),np.median(vx),np.max(vx),np.std(vx),' (min/median/max/std)')
        print('Var v_y_gas == ',np.min(vy),np.median(vy),np.max(vy),np.std(vy),' (min/median/max/std)')
        print('Var v_z_gas == ',np.min(vz),np.median(vz),np.max(vz),np.std(vz),' (min/median/max/std)')
    if('MagneticField' in P.keys()):
        bx = P['MagneticField'][:,0]; by = P['MagneticField'][:,1]; bz = P['MagneticField'][:,2]
        if(quiet==False):
            print('Var B_x_gas == ',np.min(bx),np.median(bx),np.max(bx),np.std(bx),' (min/median/max/std)')
            print('Var B_y_gas == ',np.min(by),np.median(by),np.max(by),np.std(by),' (min/median/max/std)')
            print('Var B_z_gas == ',np.min(bz),np.median(bz),np.max(bz),np.std(bz),' (min/median/max/std)')
            print('B_z_vol == ',np.sum(bz*P['Masses'][:]/P['Density'][:])/np.sum(P['Masses'][:]/P['Density'][:]))
    else:
        bx=vx; by=vy; bz=vz;
    
    zmx=np.max(zz)-np.min(zz); zzmed = np.median(zz); 
    if(zmed_set > -1.e9): zzmed=zmed_set;
    dzz=np.abs(zz-zzmed); dzz[(dzz>0.5*zmx)] = zmx-dzz[(dzz>0.5*zmx)]
    if('SmoothingLength' in P.keys()):
        ok=np.where(dzz < 0.5*P['SmoothingLength'][:])
    else:
        ok=np.where(dzz < 0.05)

    x=1.0*(xx/dmax+0.); y=1.0*(yy/dmax+0.); z=1.0*(zz/dmax+0.); 
    if((gas_val_toplot=='bvec')|(gas_val_toplot=='vvec')):
        yg,xg=np.mgrid[0:1:128j, 0:1:128j]
    else:
        yg,xg=np.mgrid[0:1:2048j, 0:1:2048j]

    v0 = P['Velocities']
    if(ptype=='PartType0'):
        u = np.log(P['Density'][:])
        u0 = P['InternalEnergy'][:]
        if(quiet==False):
            print('particle number = ',u.size)
            print('min/max/std internal energy = ',np.min(u0),np.max(u0),np.std(u0))
            print('min/max/std velocity = ',np.min(v0),np.max(v0),np.sqrt(np.std(v0[:,0])**2+np.std(v0[:,1])**2+np.std(v0[:,2])**2))
            print('min/max/std density = ',np.min(u),np.max(u),np.std(u))
    
    
    if(gas_val_toplot=='p'): 
        u = np.log( P['Density'][:] * P['InternalEnergy'][:] )
    if((gas_val_toplot=='b')|(gas_val_toplot=='bpt')|(gas_val_toplot=='bx')|(gas_val_toplot=='by')|(gas_val_toplot=='bz')):
        vx=P['MagneticField'][:,0]; vy=P['MagneticField'][:,1]; vz=P['MagneticField'][:,2];
        u=np.sqrt(vx*vx+vy*vy+vz*vz)
        if(gas_val_toplot=='bx'): u=vx
        if(gas_val_toplot=='by'): u=vy
        if(gas_val_toplot=='bz'): u=vz
    if((gas_val_toplot=='v')|(gas_val_toplot=='vpt')|(gas_val_toplot=='vx')|(gas_val_toplot=='vy')|(gas_val_toplot=='vz')):
        wt=P['Masses'][:]/np.sum(P['Masses'][:])
        vx=P['Velocities'][:,0]; vy=P['Velocities'][:,1]; vz=P['Velocities'][:,2];
        vx-=np.sum(vx*wt); vy-=np.sum(vy*wt); vz-=np.sum(vz*wt); 
        u=np.sqrt(vx*vx+vy*vy+vz*vz)
        if(gas_val_toplot=='vx'): u=vx
        if(gas_val_toplot=='vy'): u=vy
        if(gas_val_toplot=='vz'): u=vz
    if(gas_val_toplot=='bvec'):
        wt=P['Masses'][:]/np.sum(P['Masses'][:])
        bx = P['MagneticField'][:,xset]; by = P['MagneticField'][:,yset]; bz = P['MagneticField'][:,zset]
        bx-=np.sum(bx*wt); by-=np.sum(by*wt); bz-=np.sum(bz*wt); 
        u=bz;
    if(gas_val_toplot=='vvec'):
        wt=P['Masses'][:]/np.sum(P['Masses'][:])
        bx = P['Velocities'][:,xset]; by = P['Velocities'][:,yset]; bz = P['Velocities'][:,zset]
        bx-=np.sum(bx*wt); by-=np.sum(by*wt); bz-=np.sum(bz*wt); 
        u=bz;
    print('min/max_u_in_rhoplot = ',np.min(u),np.max(u))
    
    pylab.axis([0.,xmax,0.,1.])
    if(vmax==0): vmax=np.max(u)
    if(vmin==0): vmin=np.min(u)
    if(gas_val_toplot=='pt'):
        pylab.plot(x[ok],y[ok],marker=',',linestyle='',color='black',zorder=3)
        P_File.close(); return;

    dg=interpolate.griddata((x[ok],y[ok]),u[ok],(xg,yg),method='linear',fill_value=np.median(u[ok]));
    im=pylab.imshow(dg,interpolation='bicubic',vmin=vmin,vmax=vmax,cmap=cmap,extent=(0,1,0,1),zorder=1);

    if((gas_val_toplot=='vvec')|(gas_val_toplot=='bvec')):
        bxg=interpolate.griddata((x,y),bx,(xg,yg),method='linear'); 
        byg=interpolate.griddata((x,y),by,(xg,yg),method='linear');
        bg=np.sqrt(bxg*bxg+byg*byg);
        #pylab.streamplot(xg,yg,bxg,byg,color='black',linewidth=1.,arrowstyle='-',arrowsize=1,density=3.0,zorder=2)
        Q=pylab.streamplot(xg,yg,bxg,byg,color=bg,cmap='cool',density=[2.,2.],linewidth=2,arrowstyle='->',arrowsize=1.5,zorder=2)
        P_File.close(); return;






def gas_rho_image_multidim(ax, snum=3,sdir='./output',vmin=0.,vmax=0.,
        cmap='terrain',xmax=1.,xz=0,yz=0,gas_val_toplot='rho',rhocut=-2.5,save='dummy',
        zmed_set=-1.e10,quiet=False,zdir='z',zs=0.,fontsize=18.,
        center=np.zeros(0),resolution='low'):
    P_File=load_snap(sdir,snum);
    P = P_File['PartType0']
    Pc= np.array(P['Coordinates'])
    center=np.array(center); 
    if(center.size > 0): Pc-=center;
    Pc[(Pc>1.)] -= 1.;
    Pc[(Pc<0.)] += 1.;
    dmax = np.max(Pc)
    xset=0; yset=1; zset=2;
    if(xz==1): 
        xset=0; yset=2; zset=1;
    if(yz==1): 
        xset=1; yset=2; zset=0;
    xx = Pc[:,xset]; yy = Pc[:,yset]; zz = Pc[:,zset]
    vx = P['Velocities'][:,0]; vy = P['Velocities'][:,1]; vz = P['Velocities'][:,2]
    bx = P['MagneticField'][:,0]; by = P['MagneticField'][:,1]; bz = P['MagneticField'][:,2]
    if(quiet==False):
        print('Var v_x_gas == ',np.min(vx),np.median(vx),np.max(vx),np.std(vx),' (min/median/max/std)')
        print('Var v_y_gas == ',np.min(vy),np.median(vy),np.max(vy),np.std(vy),' (min/median/max/std)')
        print('Var v_z_gas == ',np.min(vz),np.median(vz),np.max(vz),np.std(vz),' (min/median/max/std)')
        print('Var B_x_gas == ',np.min(bx),np.median(bx),np.max(bx),np.std(bx),' (min/median/max/std)')
        print('Var B_y_gas == ',np.min(by),np.median(by),np.max(by),np.std(by),' (min/median/max/std)')
        print('Var B_z_gas == ',np.min(bz),np.median(bz),np.max(bz),np.std(bz),' (min/median/max/std)')
    
    zmx=np.max(zz)-np.min(zz); zzmed = np.median(zz); 
    if(zmed_set > -1.e9): zzmed=zmed_set;
    dzz=np.abs(zz-zzmed); dzz[(dzz>0.5*zmx)] = zmx-dzz[(dzz>0.5*zmx)]
    ok=np.where(dzz < 0.5*P['SmoothingLength'][:])

    x=1.0*(xx/dmax+0.); y=1.0*(yy/dmax+0.); z=1.0*(zz/dmax+0.); 
    if((gas_val_toplot=='bvec')|(gas_val_toplot=='vvec')):
        yg,xg=np.mgrid[0:1:128j, 0:1:128j]
    else:
        if(resolution=='low'): 
            yg,xg=np.mgrid[0:1:128j, 0:1:128j]
        else:
            yg,xg=np.mgrid[0:1:256j, 0:1:256j]

    u = np.log(P['Density'][:])
    labelcbar = r'$\ln{[\rho_{\rm g}/\rho_{\rm g}^{0}]}$'
    u0 = P['InternalEnergy'][:]
    v0 = P['Velocities']
    if(quiet==False):
        print('particle number = ',u.size)
        print('min/max/std internal energy = ',np.min(u0),np.max(u0),np.std(u0))
        print('min/max/std velocity = ',np.min(v0),np.max(v0),np.sqrt(np.std(v0[:,0])**2+np.std(v0[:,1])**2+np.std(v0[:,2])**2))
        print('min/max/std density = ',np.min(u),np.max(u),np.std(u))
    
    wt_vol = np.array(P['Masses'][:]/P['Density'][:])
    wt_vol[np.isfinite(wt_vol)==False]==np.median(wt_vol[np.isfinite(wt_vol)==True])
    wt_vol = wt_vol[ok]
    #wt_vol *= wt_vol;
    wt_vol /= np.sum(wt_vol)
    fillvalue = np.sum(u[ok]*wt_vol)
    if(gas_val_toplot=='p'): 
        u = np.log( P['Density'][:] * P['InternalEnergy'][:] )
        fillvalue = np.median(u[ok])
        fillvalue = np.sum(u[ok]*wt_vol)
        labelcbar = r'$\ln{[P/P_{0}]}$'
    if((gas_val_toplot=='b')|(gas_val_toplot=='bpt')):
        vx=P['MagneticField'][:,0]; vy=P['MagneticField'][:,1]; vz=P['MagneticField'][:,2];
        u=np.sqrt(vx*vx+vy*vy+vz*vz)
        fillvalue = np.median(u[ok])
        fillvalue = np.sum(u[ok]*wt_vol)
        labelcbar = r'$|{\bf B}|/\sqrt{4\pi\,P_{0}}$'
    if((gas_val_toplot=='v')|(gas_val_toplot=='vpt')|(gas_val_toplot=='vx')|(gas_val_toplot=='vy')|(gas_val_toplot=='vz')):
        wt=P['Masses'][:]/np.sum(P['Masses'][:])
        vx=P['Velocities'][:,0]; vy=P['Velocities'][:,1]; vz=P['Velocities'][:,2];
        vx-=np.sum(vx*wt); vy-=np.sum(vy*wt); vz-=np.sum(vz*wt); 
        u=np.sqrt(vx*vx+vy*vy+vz*vz)
        labelcbar = r'$|{\bf u}_{g}-\langle {\bf u}_{g} \rangle|/c_{s}^{0}$'
        if(gas_val_toplot=='vx'): 
            u=vx
            labelcbar = r'$|{\bf u}_{g,\,x}-\langle {\bf u}_{g,\,x} \rangle|/c_{s}^{0}$'
        if(gas_val_toplot=='vy'): 
            u=vy
            labelcbar = r'$|{\bf u}_{g,\,y}-\langle {\bf u}_{g,\,y} \rangle|/c_{s}^{0}$'
        if(gas_val_toplot=='vz'): 
            u=vz
            labelcbar = r'$|{\bf u}_{g,\,z}-\langle {\bf u}_{g,\,z} \rangle|/c_{s}^{0}$'
        fillvalue = np.median(u[ok])
        fillvalue = np.sum(u[ok]*wt_vol)
    if(gas_val_toplot=='bvec'):
        bx = P['MagneticField'][:,xset]; by = P['MagneticField'][:,yset]; bz = P['MagneticField'][:,zset]
        u=bz;
        fillvalue = np.median(u[ok])
        fillvalue = np.sum(u[ok]*wt_vol)
        labelcbar = r'${\bf B}_{\bot}/\sqrt{4\pi\,P_{0}}$'
    if(gas_val_toplot=='vvec'):
        wt=P['Masses'][:]/np.sum(P['Masses'][:])
        bx = P['Velocities'][:,xset]; by = P['Velocities'][:,yset]; bz = P['Velocities'][:,zset]
        bx-=np.sum(bx*wt); by-=np.sum(by*wt); bz-=np.sum(bz*wt); 
        u=bz;
        fillvalue = np.median(u[ok])
        fillvalue = np.sum(u[ok]*wt_vol)
        labelcbar = r'$({\bf u}_{g,\,\bot}-\langle {\bf u}_{g,\,\bot} \rangle)/c_{s}^{0}$'
    print('min/max_u_in_rhoplot = ',np.min(u),np.max(u))
    
    pylab.axis([0.,xmax,0.,1.])
    if(vmax==0): vmax=np.max(u)
    if(vmin==0): vmin=np.min(u)
    if(gas_val_toplot=='pt'):
        pylab.plot(x[ok],y[ok],marker=',',linestyle='',color='black')
        P_File.close(); 
        return vmin,vmax,vmax

    
    dg=interpolate.griddata((x[ok],y[ok]),u[ok],(xg,yg),method='linear',
        fill_value=fillvalue);
    dg[(dg>vmax)]=vmax; dg[(dg<vmin)]=vmin;
    cmmp=pylab.get_cmap(cmap)
    norm=matplotlib.colors.Normalize(vmin=vmin,vmax=vmax)
    colors=cmmp(norm(dg))
    if(zdir=='x'): surf=ax.plot_surface(np.zeros_like(xg)+zs,xg,yg,cstride=1,rstride=1,facecolors=colors,shade=False)
    if(zdir=='y'): surf=ax.plot_surface(xg,np.zeros_like(xg)+zs,yg,cstride=1,rstride=1,facecolors=colors,shade=False)
    if(zdir=='z'): surf=ax.plot_surface(xg,yg,np.zeros_like(xg)+zs,cstride=1,rstride=1,facecolors=colors,shade=False)
    surface_saver = norm

    if(zdir=='z'): 
        m = matplotlib.cm.ScalarMappable(cmap=cmap, norm=surf.norm)
        m.set_array(dg)
        m.set_clim(vmin=vmin,vmax=vmax)
        cb1=pylab.colorbar(m,shrink=0.7,aspect=30,fraction=0.05,pad=-0.04)
        cb1.set_label(labelcbar,fontsize=fontsize*1.15,labelpad=10.)
        cb1.ax.tick_params(labelsize=fontsize)
        cb1.ax.yaxis.set_tick_params(direction='in')

    if((gas_val_toplot=='vvec')|(gas_val_toplot=='bvec')):
        bxg=interpolate.griddata((x,y),bx,(xg,yg),method='linear'); 
        byg=interpolate.griddata((x,y),by,(xg,yg),method='linear');
        bg=np.sqrt(bxg*bxg+byg*byg);

        cmap_stream = 'cool'; cmap_final = 'magenta'
        fig_2d, ax_2d, = pylab.subplots();
        render = ax_2d.streamplot(xg,yg,bxg,byg,color=bg,cmap=cmap_stream,
            density=[2.,2.],linewidth=2,arrowstyle='->',arrowsize=1.5)
        ax_2d.clear()
        for line in render.lines.get_paths():
            new_x = line.vertices.T[0]
            new_y = line.vertices.T[1]
            z_height = np.zeros_like(new_x) + zs
            #bxg=interpolate.griddata((x,y),bx,(new_x,new_y),method='linear'); 
            #byg=interpolate.griddata((x,y),by,(new_x,new_y),method='linear');
            #bg=np.sqrt(bxg*bxg+byg*byg);
            if(zdir=='x'): 
                ax.plot(z_height,new_x,new_y,color=cmap_final,linewidth=2.)
                #ax.plot(z_height,new_x,new_y,color=bg,cmap=cmap_stream)
            if(zdir=='y'): 
                ax.plot(new_x,z_height,new_y,color=cmap_final,linewidth=2.5)
                #ax.plot(new_x,z_height,new_y,color=bg,cmap=cmap_stream)
            if(zdir=='z'): 
                ax.plot(new_x,new_y,z_height,color=cmap_final,linewidth=2.)
                #ax.plot(new_x,new_y,z_height,color=bg,cmap=cmap_stream)
        plot.close(fig_2d)


    P_File.close(); 
    return vmin,vmax,surface_saver
