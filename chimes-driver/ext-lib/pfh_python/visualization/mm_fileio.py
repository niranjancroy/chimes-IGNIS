import numpy as np
import gc
from mm_utilities import *


## here's the grunt work of loading the data we'll need 
def load_pos_vel_etc( snapdir, snapnum, filename_prefix = '', \
        h0=1, four_char=0, cosmological=0, skip_bh=0, \
        use_rundir=0, #""" make sure to change this as appropriate for the system!"""
        GAS=0, SWAP_TEMP_RHO=0, BOX_MAX=0, CENTER_BOX=[0.,0.,0.], min_stellar_age=0. ):

    ## decide which particle types to load
    have=1; ptypes=[0]; dum=np.zeros(1);
    if(GAS==0): 
        ptypes=[4]
        if (cosmological==0): ptypes=[2,3,4]

    ## check whether snapshot file exists
    fname,fname_base,fname_ext = check_if_filename_exists(snapdir,snapnum,\
        snapshot_name='snapshot',extension=extension,four_char=four_char)
    if((fname=='NULL')|(fname_ext!='.hdf5')): return dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum;
    gc.collect()

    ## load header information 
    file = h5py.File(fname,'r') # Open hdf5 snapshot file
    header_master = file["Header"] # Load header dictionary (to parse below)
    header_toparse = header_master.attrs
    numfiles = header_toparse["NumFilesPerSnapshot"]
    npartTotal = header_toparse["NumPart_Total"]
    if(np.sum(npartTotal[ptypes])<1): return dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum;
    ascale = time = header_toparse["Time"]
    hubble = header_toparse["HubbleParam"]
    if(cosmological==0): 
        ascale = 1;
    else:
        time = cosmological_time(ascale,h=hubble)
    file.close()
    r_code_to_phys = ascale / hubble; # scale to code units
    cen = np.array(CENTER_BOX)
    d_thold = 1.e10
    if(BOX_MAX != 0): d_thold = BOX_MAX / r_code_to_phys
    gc.collect()
    
    ## start parsing the actual snapshot
    xyz=np.zeros((0,3)); vxyz=xyz; q=np.zeros((0)); h_all=q; c_all=q; m_all=q; zm_all=q; id_all=np.zeros((0),dtype=int8)
    for i_file in range(numfiles):
        if (numfiles>1): fname = fname_base+'.'+str(i_file)+fname_ext  
        if(os.stat(fname).st_size>0):
            file = h5py.File(fname,'r') # Open hdf5 snapshot file
            npart = file["Header"].attrs["NumPart_ThisFile"]
            for ptype in ptypes:
                if(npart[ptype] > 1):
                    p_name_0 = 'PartType'+str(ptype)+'/'
                    xyz_all = np.array(file[p_name_0+'Coordinates'])
                    d = np.amax(np.abs(xyz_all-cen),axis=1)
                    ok = np.where(d < d_thold)[0]
                    n0 = ok.shape[0]
                    if(n0 > 0): 
                        xyz = np.concatenate([xyz,xyz_all.take(ok,axis=0)])
                        vxyz = np.concatenate([vxyz,np.array(file[p_name_0+'Velocities']).take(ok,axis=0)])
                        m_all = np.concatenate([m_all,np.array(file[p_name_0+'Masses']).take(ok)])
                        zm_all = np.concatenate([zm_all,np.array(file[p_name_0+'Metallicity']).take(ok,axis=0)[:,0]])
                        id_all = np.concatenate([id_all,np.array(file[p_name_0+'ParticleIDs']).take(ok)])
                        if(GAS==1):
                            h_all = np.concatenate([h_all,np.array(file[p_name_0+'SmoothingLength']).take(ok)])
                            if(SWAP_TEMP_RHO==1):
                                c_all = np.concatenate([c_all,np.array(file[p_name_0+'Density']).take(ok)])
                            else:
                                gas_u = np.array(file[p_name_0+'InternalEnergy']).take(ok)
                                gas_ne = np.array(file[p_name_0+'ElectronAbundance']).take(ok)
                                gas_t = 106.264 * gas_u / (1.07895 + gas_ne)
                                c_all = np.concatenate([c_all,gas_t])
                        else:
                            a0 = np.array(file[p_name_0+'StellarFormationTime']).take(ok)
                            if(cosmological==0):
                                dt = time - a0
                            else:
                                dt = time - cosmological_time(a0,h=hubble)
                            dt = np.minimum(dt , min_stellar_age)
                            c_all = np.concatenate([c_all,dt])
                    gc.collect()
                gc.collect()
            file.close()
        gc.collect()
    gc.collect()
    xyz = (xyz-cen) * ascale / hubble
    m_all /= hubble
    vxyz *= np.sqrt(ascale)
    if(GAS==1):
        h_all *= ascale / hubble
        if(SWAP_TEMP_RHO==1): c_all *= hubble*hubble / (ascale*ascale*ascale)
    else:
        h_all = load_allstars_hsml_formovie(snapdir,snapnum,cosmo=cosmological, \
            use_rundir=use_rundir,four_char=four_char,use_h0=h0,
            filename_prefix=filename_prefix,xyzset=xyz);

    gc.collect()
    return id_all, m_all, xyz[:,0], xyz[:,1], xyz[:,2], vxyz[:,0], vxyz[:,1], vxyz[:,2], 1.25*h_all, c_all, zm_all;


    
def load_allstars_hsml_formovie(snapdir,snapnum,cosmo=0,use_rundir=0,
        four_char=0,use_h0=1,filename_prefix='',xyzset=np.array([0.,0.,0.])):

    rootdir='/work/01799/phopkins/stellar_hsml/'
    exts=snap_ext(snapnum,four_char=four_char);
    s0=snapdir.split("/"); ss=s0[len(s0)-1]; 
    if(len(ss)==0): ss=s0[len(s0)-2];
    if(filename_prefix==''):
        hsmlfile_r=rootdir+ss+'_allstars_hsml_'+exts
    else:
        hsmlfile_r=rootdir+filename_prefix+'_allstars_hsml_'+exts
    if (use_rundir==1): hsmlfile_r=snapdir+'/allstars_hsml_'+exts
    hsmlfile=hsmlfile_r+'.dat' ## check if binary file exists

    if os.path.exists(hsmlfile): ## it exists! 
        lut=open(hsmlfile,'r'); 
        int_in=array.array('i'); int_in.fromfile(lut,1); nstars=int_in[0];
        h_in=array.array('f'); h_in.fromfile(lut,nstars);         
        lut.close();
        gc.collect()
        return np.array(np.copy(h_in));
    else: ## no pre-computed file, need to do it ourselves
        h = get_particle_hsml(xyzset[:,0],xyzset[:,1],xyzset[:,2],DesNgb=62);
        ## great now we've got the stars, lets write this to a file for next time
        nstars = checklen(h);
        print(nstars)
        if (nstars>1):
            lut=open(hsmlfile,'wb');
            lut.write(struct.pack('i',nstars));
            lut.write(struct.pack('f'*len(h),*h));
            lut.close();
            gc.collect()
            return h;
    gc.collect()
    return 0;



def get_particle_hsml( x, y, z, DesNgb=32, Hmax=0.):
    x=fcor(x); y=fcor(y); z=fcor(z); N=checklen(x); 
    ok=(ok_scan(x) & ok_scan(y) & ok_scan(z)); x=x[ok]; y=y[ok]; z=z[ok];
    if(Hmax==0.):
        dx=np.max(x)-np.min(x); dy=np.max(y)-np.min(y); dz=np.max(z)-np.min(z); ddx=np.max([dx,dy,dz]); 
        Hmax=5.*ddx*(np.float(N)**(-1./3.)); ## mean inter-particle spacing

    ## load the routine we need
    exec_call=util.return_python_routines_cdir()+'/StellarHsml/starhsml.so'
    h_routine=ctypes.cdll[exec_call];

    h_out_cast=ctypes.c_float*N; H_OUT=h_out_cast();
    ## main call to the hsml-finding routine
    h_routine.stellarhsml( ctypes.c_int(N), \
        vfloat(x), vfloat(y), vfloat(z), ctypes.c_int(DesNgb), \
        ctypes.c_float(Hmax), ctypes.byref(H_OUT) )

    ## now put the output arrays into a useful format 
    h = np.ctypeslib.as_array(np.copy(H_OUT));
    gc.collect()
    return h;





##
## handy function to load everything we need from snapshots, in particular to 
##   concatenate the various star types if we're using non-cosmological snapshots
##
def load_snapshot_brick(ptypes, snapdir, snapnum, h0=1, cosmological=0, \
        skip_bh=0, do_xray=0, four_char=0, use_rundir=0, full_gas_set=0 , filename_prefix=''):
    have=0; have_h_stars=0;
    ppp_head=gadget.readsnap(snapdir,snapnum,0,h0=h0,cosmological=cosmological,skip_bh=skip_bh,header_only=1)
    time=ppp_head['time']
    if (full_gas_set==1):
        ## here it will just return the entire gas 'brick'
        ppp=gadget.readsnap(snapdir,snapnum,0,h0=h0,cosmological=cosmological,skip_bh=skip_bh);
        m=ppp['m']; p=ppp['p']; x=p[:,0]; y=p[:,1]; z=p[:,2]; zm=ppp['z']; 
        if(checklen(zm.shape)>1): zm=zm[:,0]; lx=get_gas_xray_luminosity(ppp) / (3.9e33);
        sfr0=ppp['SFR']; sfr=1.*sfr0; 
        ok=sfr0[(sfr0 > 0.)]
        sfr = 1.e-10 * 0.25*(1./(1.+(ppp['rho'])**(-0.5)))
        if(ok.size > 0):
            p0=np.min(ppp['rho'][(sfr0 > 0.)]); 
            lo=(sfr0 <= 0.); sfr[lo]=0.25*(1./(1.+(ppp['rho'][lo]/p0)**(-0.5)))*np.min(sfr0[(sfr0>0.)]);
        gc.collect()
        return ppp['u'], ppp['rho'], ppp['h'], ppp['nh'], ppp['ne'], sfr, lx, zm, m, x, y, z, time;
        
    for ptype in ptypes:
        ppp=gadget.readsnap(snapdir,snapnum,ptype,h0=h0,cosmological=cosmological,skip_bh=skip_bh);
        if(ppp['k']==1):
            n=checklen(ppp['m']);
            if(n>1):
                m=ppp['m']; p=ppp['p']; x=p[:,0]; y=p[:,1]; z=p[:,2]; 
                if (ptype==0): 
                    cc=gadget.gas_temperature(ppp['u'],ppp['ne']);
                    if (do_xray==1): m=get_gas_xray_luminosity(ppp) / (3.9e33); ## bremstrahhlung (in solar)
                    h=ppp['h']; zm=ppp['z']; 
                    if(checklen(zm.shape)>1): zm=zm[:,0]
                if (ptype==4):
                    cc=gadget.get_stellar_ages(ppp,ppp_head,cosmological=cosmological);
                    zm=ppp['z']; 
                    if(checklen(zm.shape)>1): zm=zm[:,0]
                if (ptype==2): ## need to assign ages and metallicities
                    cc=np.random.rand(n)*(ppp_head['time']+4.0);
                    zm=(np.random.rand(n)*(1.0-0.1)+0.1) * 0.02;
                if (ptype==3): ## need to assign ages and metallicities
                    cc=np.random.rand(n)*(ppp_head['time']+12.0);
                    zm=(np.random.rand(n)*(0.3-0.03)+0.03) * 0.02;
                hstars_should_concat=0;
                if((ptype>0) & (have_h_stars==0)):
                    h=starhsml.load_allstars_hsml(snapdir,snapnum,cosmo=cosmological, \
                        use_rundir=use_rundir,four_char=four_char,use_h0=h0,filename_prefix=filename_prefix);
                    have_h_stars=1;
                    hstars_should_concat=1;
                if (have==1):
                    m_all=np.concatenate((m_all,m)); x_all=np.concatenate((x_all,x)); 
                    y_all=np.concatenate((y_all,y)); z_all=np.concatenate((z_all,z)); 
                    c_all=np.concatenate((c_all,cc)); zm_all=np.concatenate((zm_all,zm)); 
                    if(hstars_should_concat==1): h_all=np.concatenate((h_all,h)); 
                else:
                    m_all=m; x_all=x; y_all=y; z_all=z; zm_all=zm; c_all=cc; h_all=h; 
                    have=1;
    if (have==1):
        gc.collect()
        return m_all, x_all, y_all, z_all, c_all, h_all, zm_all, time;
    gc.collect()
    return 0,0,0,0,0,0,0,0;




