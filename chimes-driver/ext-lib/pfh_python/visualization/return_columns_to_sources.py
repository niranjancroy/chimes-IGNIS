import numpy as np
import math
import ctypes
import pfh_utils as util
import gadget

def checklen(x):
    return len(np.array(x,ndmin=1));

def int_round(x):
    return np.int(np.round(x));

def ok_scan(input,xmax=1.0e10,pos=0):
    if (pos==1):
        return (np.isnan(input)==False) & (np.fabs(input)<=xmax) & (input > 0.);
    else:
        return (np.isnan(input)==False) & (np.fabs(input)<=xmax);

def fcor(x):
    return np.array(x,dtype='f',ndmin=1)
def vfloat(x):
    return x.ctypes.data_as(ctypes.POINTER(ctypes.c_float));



def calculate_zoom_center(sdir,snum,cen=[0.,0.,0.],clip_size=2.e10):
    rgrid=np.array([1.0e10,1000.,700.,500.,300.,200.,100.,70.,50.,30.,20.,10.,5.,2.5,1.]);
    rgrid=rgrid[rgrid <= clip_size];
    Ps=gadget.readsnap(sdir,snum,4,cosmological=1);
    n_new=Ps['m'].shape[0];
    if (n_new > 1):
        pos=Ps['p']; x0s=pos[:,0]; y0s=pos[:,1]; z0s=pos[:,2];
    Pg=gadget.readsnap(sdir,snum,0,cosmological=1);
    rho=np.array(Pg['rho'])*407.5;
    if (rho.shape[0] > 0):
        pos=Pg['p']; x0g=pos[:,0]; y0g=pos[:,1]; z0g=pos[:,2];
    rho_cut=1.0e-5;
    cen=np.array(cen);

    for i_rcut in range(len(rgrid)):
        for j_looper in range(5):
            if (n_new > 1000):
                x=x0s; y=y0s; z=z0s;
            else:
                ok=(rho > rho_cut);
                x=x0g[ok]; y=y0g[ok]; z=z0g[ok];
            x=x-cen[0]; y=y-cen[1]; z=z-cen[2];
            r = np.sqrt(x*x + y*y + z*z);
            ok = (r < rgrid[i_rcut]);
            if (len(r[ok]) > 1000):
                x=x[ok]; y=y[ok]; z=z[ok];
                if (i_rcut <= len(rgrid)-5):
                    cen+=np.array([np.median(x),np.median(y),np.median(z)]);
                else:
                    cen+=np.array([np.mean(x),np.mean(y),np.mean(z)]);
            else:
                if (len(r[ok]) > 200):
                    x=x[ok]; y=y[ok]; z=z[ok];
                    cen+=np.array([np.mean(x),np.mean(y),np.mean(z)]);
                    
    return cen;
    

def compare_los_columns(filename='columns'):
    import numpy as np
    import h5py as h5py
    import os.path
    import math
    import matplotlib
    import cosmology as my_cosmo
    matplotlib.use('Agg') ## this calls matplotlib without a X-windows GUI
    import matplotlib.pylab as pylab
    import matplotlib.pyplot as plot
    import scipy.fftpack as fft
    import scipy.misc
    import scipy.interpolate as interpolate
    import scipy.optimize as optimize
    import scipy.special as scifun
    import gadget
    import pfh_utils as util
    import halo_finder as halofind
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    from scipy import ndimage

    pylab.close('all'); plot.figure(1,figsize=(5.1,3.4)); f00=16.; f00t=f00;
    matplotlib.rcParams.update({'font.size':f00})

    F=h5py.File(filename+'.dat')
    sdir=F['Snapshot_Directory'].value; snum=F['Snapshot_Number'].value; cen=F['Snapshot_Center'].value; xlen=F['Side_Length_Box'].value;
    nlos=F['N_Line_Of_Sight'].value; phi_vec=F['Phi_List'].value; cos_theta_vec=F['Cos_Theta_List'].value; agemax=F['Max_Age_Used'].value;
    los_nh=F['LOS_nh_column'].value; los_z=F['LOS_Z_weighted'].value;
    
    print(F['LOS_nh_column'].value)
    print(F['LOS_Z_weighted'].value)
    
    PPPh=gadget.readsnap(sdir,snum,4,cosmological=1,header_only=1);
    PPPs=gadget.readsnap(sdir,snum,4,cosmological=1);
    age = gadget.get_stellar_ages(PPPs,PPPh,cosmological=1);
    source_pos=PPPs['p']; source_pos[:,0]-=cen[0]; source_pos[:,1]-=cen[1]; source_pos[:,2]-=cen[2]; 
    r_s=np.sqrt(source_pos[:,0]**2+source_pos[:,1]**2+source_pos[:,2]**2)
    ok_s=(age<agemax)&(r_s<xlen);

    F0=h5py.File(sdir+'/snapshot_'+str(snum)+'.hdf5','r')
    Fs=F0['PartType4']
    Z_s = Fs['Metallicity'][:,0]
    NH_noH = Fs['Metallicity'][:,1]
    NH_wH = Fs['Metallicity'][:,2]
    L_s = Fs['Metallicity'][:,3]
    
    code_unit_to_NH = 8.741e23
    
    okk = (age[ok_s] < 0.25)
    
    NH_code_Z = (1./0.7) * NH_wH[ok_s][okk] * (Z_s[ok_s][okk]/0.02)
    NH_post_Z = los_nh[okk,:] * (los_z[okk,:]/0.02) + 1.e14
    
    kappa = 1.67e-24 * 180.
    q = NH_post_Z * kappa
    t = np.exp(-q)
    x = np.log10(NH_post_Z)
    tx = t*x
    ts=np.sum(t,axis=1); txs=np.sum(tx,axis=1); tmean=txs/ts
    
    Lwt = L_s[ok_s][okk]
    
    fabs = np.sum(1.-t,axis=1)/(1.*nlos)
    print('MINMAX == ',np.max(fabs),np.min(fabs))
    fabs_code = (1.-np.exp(-NH_code_Z*code_unit_to_NH*kappa))
    print('MINMAX == ',np.max(fabs_code),np.min(fabs_code))
    
    x = np.arange(19.,24.,0.1)
    pylab.plot(x,x,linestyle='-',linewidth=2.,color='white')
    for i in range(nlos):
        x = np.log10(NH_code_Z * code_unit_to_NH)
        y = np.log10(NH_post_Z[:,i])
        #pylab.plot(x,y,marker=',',linestyle='',rasterized=True,color='blue')

        #x = np.log10(NH_wH[ok_s][okk] * code_unit_to_NH)
        #y = np.log10(los_nh[okk,i] + 1.e14)
        #pylab.plot(x,y,marker=',',linestyle='',rasterized=True,color='blue')

        pylab.hexbin(x,y,gridsize=(50,50),cmap='nipy_spectral',bins='log',
            extent=(18.,24.,18.,24.))

    pylab.axis([19.,24.,19.,24.])
    F.close()
    F0.close()

    pylab.close('all'); plot.figure(1,figsize=(5.1,3.4)); f00=16.; f00t=f00;
    matplotlib.rcParams.update({'font.size':f00})
    
    pylab.plot(np.log10(NH_code_Z * code_unit_to_NH), tmean,
        color='blue',marker=',',linestyle='')
    x = np.arange(19.,24.,0.1)
    pylab.plot(x,x,linestyle='-',linewidth=2.,color='black')
    

    if(2==0):
        pylab.close('all'); plot.figure(1,figsize=(5.1,3.4)); f00=16.; f00t=f00;
        matplotlib.rcParams.update({'font.size':f00})
        bins = 100
        wtf = 0.*Lwt+1.
        pylab.yscale('log')
    
        pylab.hist(np.log10(fabs),bins=bins,weights=wtf,normed=True,histtype='step',
            color='black',cumulative=0)
        pylab.hist(np.log10(fabs_code),bins=bins,weights=wtf,normed=True,histtype='step',
            color='red',cumulative=0)

    if(2==0):
        pylab.close('all'); plot.figure(1,figsize=(5.1,3.4)); f00=16.; f00t=f00;
        matplotlib.rcParams.update({'font.size':f00})
        bins = 50
        wtf = 1.*Lwt+0.
        pylab.yscale('linear'); pylab.axis([-1.5,1.5,0.,14.2])
        pylab.yscale('log'); pylab.axis([-1.5,1.5,0.01,14.2])
        pylab.xlabel(r'$\log{(\ P_{\rm rad}^{\rm LEBRON}\ / P_{\rm rad}^{\rm Full-RT}\ )}$')
        pylab.ylabel(r'${\rm Fraction}$')

        x = 0.9*(np.log10(fabs_code)-np.log10(fabs))
        w = wtf*fabs_code
        pylab.hist(x,bins=bins,weights=w,normed=True,histtype='step',color='black')
        #pylab.hist(x,bins=bins,weights=w,normed=True,histtype='step',cumulative=-1,color='red')



    
def dump_los_columns_vs_angle(xlen=10.,n_los=25,filename='columns',agemax=0.25):
    import gadget
    import h5py
    sdir='/Users/phopkins/Downloads/gal_lookie/m12i_ref12_560_StellarNHSaved'; snum=561; cen=[ 39802.78755785,41916.7327608,43891.70273863]
    #sdir='/Users/phopkins/Downloads/gal_lookie/m12i_ref13'; snum=600; cen=[ 41875.75676213,  44122.38220725,  46257.48356735]

    PPPh=gadget.readsnap(sdir,snum,4,cosmological=1,header_only=1);
    PPPs=gadget.readsnap(sdir,snum,4,cosmological=1);
    PPP=gadget.readsnap(sdir,snum,0,cosmological=1);
    age = gadget.get_stellar_ages(PPPs,PPPh,cosmological=1);
    source_pos=PPPs['p']; gas_pos=PPP['p'];
    source_pos[:,0]-=cen[0]; gas_pos[:,0]-=cen[0]; source_pos[:,1]-=cen[1]; gas_pos[:,1]-=cen[1]; source_pos[:,2]-=cen[2]; gas_pos[:,2]-=cen[2]; 
    r_s=np.sqrt(source_pos[:,0]**2+source_pos[:,1]**2+source_pos[:,2]**2)
    r_g=np.sqrt(gas_pos[:,0]**2+gas_pos[:,1]**2+gas_pos[:,2]**2)
    ok_s=(age<agemax)&(r_s<xlen); ok_g=(r_g<xlen); 
    source_pos=source_pos[ok_s,:]; gas_pos=gas_pos[ok_g,:]
    
    cos_theta_vec=-1.+2.*np.random.rand(n_los)
    phi_vec=2.*np.pi*np.random.rand(n_los)
    source_pos_0=source_pos; gas_pos_0=gas_pos;
    n_sources = r_s[ok_s].size
    n_angles = cos_theta_vec.size
    b_arr = np.zeros((n_sources,n_angles))
    los_nh_all = np.zeros((n_sources,n_angles)); 
    los_z_all = np.zeros((n_sources,n_angles));
    print('n_sources=',n_sources,' n_angles=',n_angles)
    for cos_theta,phi,i_ang in zip(cos_theta_vec,phi_vec,range(n_angles)):
        print('cos_theta=',cos_theta,' phi=',phi)
        ct=cos_theta; st=np.sqrt(1.-ct*ct);
        cp=np.cos(phi); sp=np.sin(phi);
        p=source_pos_0; x=p[:,0]; y=p[:,1]; z=p[:,2]; xp=x; yp=ct*y+st*z; zp=-st*y+ct*z; xq=cp*xp+sp*yp; yq=-sp*xp+cp*yp; zq=zp; p[:,0]=xq; p[:,1]=yq; p[:,2]=zq; 
        source_pos=1.*p;
        p=gas_pos_0; x=p[:,0]; y=p[:,1]; z=p[:,2]; xp=x; yp=ct*y+st*z; zp=-st*y+ct*z; xq=cp*xp+sp*yp; yq=-sp*xp+cp*yp; zq=zp; p[:,0]=xq; p[:,1]=yq; p[:,2]=zq; 
        gas_pos=1.*p;
        
        print(source_pos.shape, gas_pos.shape)

        los_nh,los_nh_hot,los_z = \
          return_columns_to_sources(source_pos,\
          gas_pos,PPP['u'][ok_g],PPP['rho'][ok_g],PPP['h'][ok_g],PPP['ne'][ok_g],PPP['nh'][ok_g],PPP['z'][ok_g],PPP['m'][ok_g],\
          xrange=[-xlen,xlen],yrange=[-xlen,xlen],zrange=[-xlen,xlen])
        print(los_nh)
        print(los_z)
        print(' ')
        print(los_nh.shape)
        print(los_nh_all.shape)
        print(i_ang)
        print(los_nh[:])
        print(los_z[:])
        los_nh_all[:,i_ang] = 1.0 * los_nh[:]
        los_z_all[:,i_ang] = 1.0 * los_z[:]
        
    print(los_nh_all)
    outfi = h5py.File(filename+'.dat','w')
    dset_sdir = outfi.create_dataset('Snapshot_Directory',data=sdir)
    dset_snum = outfi.create_dataset('Snapshot_Number',data=snum)
    dset_cen = outfi.create_dataset('Snapshot_Center',data=cen)
    dset_xlen = outfi.create_dataset('Side_Length_Box',data=xlen)
    dset_mage = outfi.create_dataset('Max_Age_Used',data=agemax)
    dset_nlos = outfi.create_dataset('N_Line_Of_Sight',data=n_los)
    dset_ctheta = outfi.create_dataset('Cos_Theta_List',data=cos_theta_vec)
    dset_phi = outfi.create_dataset('Phi_List',data=phi_vec)
    dset_los_nh = outfi.create_dataset('LOS_nh_column',data=los_nh_all)
    dset_los_z = outfi.create_dataset('LOS_Z_weighted',data=los_z_all)
    outfi.close()



##
## return: los_NH_allgas, los_NH_hotphase, los_gas_metallicity 
##
def return_columns_to_sources( source_pos, gas_pos, \
    gas_u, gas_rho, gas_hsml, gas_numh, gas_nume, gas_metallicity, gas_mass, \
    xrange=0, yrange=0, zrange=0, \
    MIN_CELL_SIZE=0.01, OUTER_RANGE_OF_INT=1200., \
    TRIM_PARTICLES=1 ):
    
    ## check the ordering of the position matrices:
    if ((checklen(gas_pos[0,:])==3) & (checklen(gas_pos[:,0]) !=3)): gas_pos=np.transpose(gas_pos);
    if ((checklen(source_pos[0,:])==3) & (checklen(source_pos[:,0]) !=3)): source_pos=np.transpose(source_pos);
    ## and that metallicities are a vector, not a matrix
    if (len(gas_metallicity.shape)>1): gas_metallicity=gas_metallicity[:,0]

    if ((checklen(gas_pos[:,0]) != 3) | (checklen(gas_pos[0,:]) <= 1)):
        print('ERROR WILL OCCUR :: need pos to be (3,N)')

    x=source_pos[0,:] ; y=source_pos[1,:] ; z=source_pos[2,:]
    if(checklen(xrange)<=1): xrange=[np.min(x),np.max(x)];
    if(checklen(yrange)<=1): yrange=[np.min(y),np.max(y)];
    xr=xrange; yr=yrange;
    if(checklen(zrange)<=1):
        zrr=np.sqrt((xr[1]-xr[0])**2.+(yr[1]-yr[0])**2.)/np.sqrt(2.);
        zmin=np.median(z)-zrr; zmax=np.median(z)+zrr;
        if (np.min(z) > zmin): zmin=np.min(z);
        zrange=[zmin,zmax]; print('z_range (calc) == ',zrange)
    zr=zrange;
    x00=0.5*(xr[1]+xr[0]); y00=0.5*(yr[1]+yr[0]); z00=0.5*(zr[1]+zr[0]); 
    tolfac = 1.0e10;
    if (TRIM_PARTICLES==1):
        tolfac = 0.05; 
        #tolfac = -0.01;
        ## trim down the incoming list to only whats in the range plotted 
        ##   (saves a ton of time and memory overflow crashes)

    dx=(0.5+tolfac)*(xr[1]-xr[0]); dy=(0.5+tolfac)*(yr[1]-yr[0]); dz=(0.5+tolfac)*(zr[1]-zr[0]);
    ok_sources=ok_scan(x-x00,xmax=dx) & ok_scan(y-y00,xmax=dy) & ok_scan(z-z00,xmax=dz);
    x=gas_pos[0,:] ; y=gas_pos[1,:] ; z=gas_pos[2,:]
    gw=gas_rho ; gh=gas_hsml ; gz=gas_metallicity ; gm=gas_mass
    ok_gas=ok_scan(x-x00,xmax=dx) & ok_scan(y-y00,xmax=dy) & ok_scan(z-z00,xmax=dz) & \
        ok_scan(gw,pos=1) & ok_scan(gh,pos=1) & ok_scan(gz,pos=1) & ok_scan(gm,pos=1,xmax=1.0e40);

    Ngas = checklen(gas_mass[ok_gas]);
    Nstars = checklen(source_pos[0,ok_sources]);
    if (Nstars<=1) or (Ngas<=1):
        print(' UH-OH: EXPECT ERROR NOW, there are no valid source/gas particles to send!')
        print('Ngas=',Ngas,'Nstars=',Nstars,'dx=',dx,'dy=',dy,'dz=',dz,'x00=',x00,'y00=',y00,'z00=',z00)
        return -1,-1,-1;

    dzmax=np.max(gas_pos[2,ok_gas])-z00; 
    if(dzmax<OUTER_RANGE_OF_INT): OUTER_RANGE_OF_INT=dzmax;
    print('PASSING: N_gas=',Ngas,'N_sources=',Nstars,'MaxDist=',OUTER_RANGE_OF_INT,'MinCell=',MIN_CELL_SIZE)
    Nbh=0; theta=1.0e-4; phi=1.0e-4;
  
    ## load the routine we need
    exec_call=util.return_python_routines_cdir()+'/LOS_column_singlePOV/getnh.so'
    NH_routine=ctypes.cdll[exec_call];
    ## cast the variables to store the results
    nh_out_cast=ctypes.c_float*Nstars; 
    los_NH_out=nh_out_cast(); los_NH_hot_out=nh_out_cast(); los_Z_out=nh_out_cast();

    ## ok this is a bit arcane but the routine will read appropriately this block order
    Coord = np.zeros((Ngas+Nstars,10),dtype='f');
    Coord[0:Ngas,0] = gas_pos[0,ok_gas]-x00;
    Coord[0:Ngas,1] = gas_pos[1,ok_gas]-y00;
    Coord[0:Ngas,2] = gas_pos[2,ok_gas]-z00;
    Coord[0:Ngas,3] = gas_u[ok_gas]
    Coord[0:Ngas,4] = gas_rho[ok_gas]
    Coord[0:Ngas,5] = gas_hsml[ok_gas]
    Coord[0:Ngas,6] = gas_numh[ok_gas]
    Coord[0:Ngas,7] = gas_nume[ok_gas]
    Coord[0:Ngas,8] = gas_metallicity[ok_gas]
    Coord[0:Ngas,9] = gas_mass[ok_gas]
    Coord[Ngas:Nstars+Ngas,0] = source_pos[0,ok_sources]-x00;
    Coord[Ngas:Nstars+Ngas,1] = source_pos[1,ok_sources]-y00;
    Coord[Ngas:Nstars+Ngas,2] = source_pos[2,ok_sources]-z00;
    Coord=np.copy(np.transpose(Coord));

    ## main call to the NH-calculation routine
    NH_routine.getnh( ctypes.c_int(Ngas), ctypes.c_int(Nstars), ctypes.c_int(Nbh), \
        ctypes.c_float(theta), ctypes.c_float(phi), \
        vfloat(Coord), \
        ctypes.byref(los_NH_out),  ctypes.byref(los_NH_hot_out),  ctypes.byref(los_Z_out), \
        ctypes.c_float(OUTER_RANGE_OF_INT), ctypes.c_float(MIN_CELL_SIZE) );
    ## now put the output arrays into a useful format 
    print(type(los_NH_out), los_NH_out)
    los_NH = np.ctypeslib.as_array(np.copy(los_NH_out));
    los_NH_hot = np.ctypeslib.as_array(np.copy(los_NH_hot_out));
    los_Z = np.ctypeslib.as_array(np.copy(los_Z_out));

    # trap for really low NH value and zero metallicity (make it small instead)
    low_NH = 1.0e10;
    los_NH[los_NH<low_NH]=low_NH; los_NH_hot[los_NH_hot<low_NH]=low_NH;
    los_Z[los_Z<=1.0e-5]=1.0e-5;

    ## assign strong attenuation to all 'off-grid' sources, then fill in calc. vals
    Nstarstot=checklen(source_pos[0,:]);
    los_NH_allgas=np.zeros(Nstarstot,dtype='f')+1.0e23;
    los_NH_hotgas=np.zeros(Nstarstot,dtype='f')+1.0e23;
    los_gas_metallicity=np.zeros(Nstarstot,dtype='f')+0.02;
    nok=checklen(los_NH_allgas[ok_sources])
    los_NH_allgas[ok_sources]=fcor(los_NH[0:Nstars]);
    los_NH_hotgas[ok_sources]=fcor(los_NH_hot[0:Nstars]);
    los_gas_metallicity[ok_sources]=fcor(los_Z[0:Nstars]);

    return los_NH_allgas, los_NH_hotgas, los_gas_metallicity;
