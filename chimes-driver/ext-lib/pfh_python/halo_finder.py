import numpy as np
import math
import ctypes
import os.path
import pfh_utils as util
import gadget
import h5py
import sys

def testg():
    smaster='/Users/phopkins/Documents/work/plots/zooms/'
    sdir='m14_tst'
    snum=275
    minmass=0.3e11
    
    smaster='/Users/phopkins/Downloads/m12_zoom/'
    sdir='z5_box'
    snum=4
    minmass=8.e9

    smaster='/Users/phopkins/Downloads/m12_zoom/'
    sdir='z5_box11'
    snum=4
    minmass=1.e10
    compile_halo_properties(smaster+sdir, snum, min_halomass_res=minmass/1.0e10)



def compile_halo_properties(snapdir, snapnum, \
        cosmological=1,h0=1,skip_bh=1,four_char=0, \
        min_halomass_res=1.0e-20):

    dm_only = 1
    clip_for_speed = 1 ## enable several lines below to limit halos scanned 
    if (dm_only==1): clip_for_speed=0


    ## first load the halo ID set: 
    halo_ids_dm = load_halo_group_ids(snapdir,snapnum, \
            cosmological=cosmological,h0=h0,skip_bh=skip_bh,four_char=four_char)
    
    ## load all the data sets for the snapshot (so we don't have to re-do this in-loop):
    ## first collate all the DM particles:
    have_key = 0
    j_check=0;
    while (have_key==0):
        P_header = gadget.readsnap(snapdir,snapnum,j_check,header_only=1,\
            h0=h0,cosmological=cosmological,skip_bh=skip_bh,four_char=four_char)
        if (P_header['k'] != -1): have_key=1;
        j_check += 1;
        if (j_check > 5): have_key=1;
    have_key = 0
    for ptype in [1,2,3,5]:
        P=gadget.readsnap(snapdir,snapnum,ptype,\
            h0=h0,cosmological=cosmological,skip_bh=skip_bh,four_char=four_char)
        if (P['k']==1): 
            if (ptype<5) or (P['m'].shape[0]>1000):
                if (have_key==0):
                    pos=P['p']; m=P['m']; type=np.zeros(P['m'].size)+np.float(ptype)
                    have_key = 1
                else:
                    pos=np.concatenate((pos,P['p'])); 
                    m=np.concatenate((m,P['m']));
                    type=np.concatenate((type,np.zeros(P['m'].size)+np.float(ptype)));
    m_dm=m; x_dm=pos[:,0]; y_dm=pos[:,1]; z_dm=pos[:,2]; type_dm=type;
    m_dm_hires_particle = np.median(m_dm[type_dm==1])
    pmin=np.min(pos); pmax=np.max(pos); box_size=pmax-pmin;
    ## and now just load the gas & star particle types
    PPP_gas = gadget.readsnap(snapdir,snapnum,0,\
            h0=h0,cosmological=cosmological,skip_bh=skip_bh,four_char=four_char)
    PPP_stars = gadget.readsnap(snapdir,snapnum,4,\
            h0=h0,cosmological=cosmological,skip_bh=skip_bh,four_char=four_char)


    ## ok now loop through and flag everything we care about in the files:
    min_number_of_particles_in_halo = 100 
    min_mass_of_halo = 2. * np.float(min_number_of_particles_in_halo) * m_dm_hires_particle
    if (clip_for_speed==1):
        if (min_halomass_res>min_mass_of_halo): min_mass_of_halo=min_halomass_res
    N_stars_in_halo_min = 10
    if (dm_only==1):
        N_stars_total = 0
    else:
        N_stars_total = PPP_stars['m'].size
    N_stars_cumulative_processed = 0

    ## prepare the dictionaries::
    hhh_dm = np.empty(0,dtype=object)
    hhh_star = np.empty(0,dtype=object)
    hhh_gas = np.empty(0,dtype=object)

    halo_id_max = max(halo_ids_dm)
    print('N_HALOS == ',halo_id_max)
    print('.. beginning loop now...')
    for halo_id in range(halo_id_max+1):
        in_halo = ( halo_ids_dm == halo_id )
        m = m_dm[in_halo]
    
        ## halos are already sorted by total particle number, so once we get to those 
        ##   which are 'too small', we can terminate this loop
        N_dm_halo = m.size
        if (N_dm_halo < min_number_of_particles_in_halo): break
        m_dm_halo = np.sum(m)
        if (m_dm_halo < min_mass_of_halo): break
        
        if (clip_for_speed==1):        
            ## also stop once the halos processed account for all the stars in the simulation, 
            if (N_stars_cumulative_processed >= N_stars_total): break
    
        ## now get the halo position
        x = x_dm[in_halo]; y = y_dm[in_halo]; z = z_dm[in_halo]; wt = m/m_dm_halo;
        cen = get_halo_com(x, y, z, wt, box_size)
    
        dx=periodic_dx(x,cen[0],box_size)
        dy=periodic_dx(y,cen[1],box_size)
        dz=periodic_dx(z,cen[2],box_size)
        r2 = dx*dx + dy*dy + dz*dz
        rmax = np.sqrt(np.max(r2))
        rrms = np.sqrt(np.sum(m*r2)/m_dm_halo)
        rmed = np.sqrt(np.median(r2))
        ## not the same as half-mass, if multi-mass particles
    
        ## also get the amount of contamination 
        t = type_dm[in_halo]
        contamination = 1. - np.float(t[t==1].size) / np.float(N_dm_halo)
    
        ## now consider stars and/or gas, recalling that each unit here is a 'central galaxy' ::
        r_cut_halo = 2.0*rrms
        if(clip_for_speed==1): r_cut_halo = 1.4*rrms
        time = P_header['time']
        
        ggg_dm = {'id':halo_id,'n':N_dm_halo,'m':m_dm_halo,'pos':cen,'r_max':rmax,\
            'r_med':rmed,'r_rms':rrms,'contamination':contamination,'time':time}
        print('(( ',halo_id,' ))')
        print(ggg_dm)
        hhh_dm = np.append( hhh_dm, ggg_dm )

        if (dm_only==0):
            ggg_star = get_galaxy_properties( PPP_stars, P_header, 4, cen, r_cut_halo, box_size=box_size )
            hhh_star = np.append( hhh_star, ggg_star )
            N_stars_cumulative_processed += ggg_star['n']
            print(ggg_star)

            ## optionally clip here (skip evaluation of gas) if no stars in the halo:
            if(clip_for_speed==1): 
                if(ggg_star['n']<=0): r_cut_halo*= 0.0
            ggg_gas = get_galaxy_properties( PPP_gas, P_header, 0, cen, r_cut_halo, \
                box_size=box_size, r_eff_cut=2.*ggg_star['r_half'] )
            hhh_gas = np.append( hhh_gas, ggg_gas )
            print(ggg_gas)
        
        print(' ')
        print(' ')
        sys.stdout.flush()


    ## ok done with major loop over halos
    N_halos_in_final = hhh_dm.size
    print('... done with census: N_halos_final = ',N_halos_in_final)
    
    ## open the file to record all saved properties (so they can be easily loaded later)
    fname = get_halo_id_filename(snapdir,snapnum,full_halo_results=1,four_char=four_char)
    print('writing to file : ',fname)
    outfi = h5py.File(fname,'w')
    
    ## organize keywords :: 
    type_vec=[0,1,2]
    if (dm_only==1): type_vec=[0]
    for type_loop in type_vec:
        if (type_loop==0):
            dic = hhh_dm
            prefix = 'halo_'
        if (type_loop==1):
            dic = hhh_gas
            prefix = 'gas_'
        if (type_loop==2):
            dic = hhh_star
            prefix = 'star_'
        for key_var in dic[0].keys():
            ## fill in an array with this keyname
            d00 = np.empty(0)
            for i in range(N_halos_in_final): 
                d0_tmp = dic[i]
                d00 = np.append(d00,d0_tmp[key_var])
            d0 = np.array(d0_tmp[key_var])
            if (d0.size>1): d00=d00.reshape(-1,d0.size)
            dataset = outfi.create_dataset(prefix+key_var,data=d00)
    outfi.close() ## done writing
    print('done!')

    ## (to read) :::
    #infiname = get_halo_id_filename(snapdir,snapnum,full_halo_results=1,four_char=four_char)
    #infi=h5py.File(infiname,'r')
    ## now any keyword can be called by its name
    #sigma = np.array(infi["star_v_disp"])   # etc
    ## use 'list(infi)' to see available keys
    #infi.close()
    
    return 1



def get_halo_com(x,y,z,wt,box):
    ## problem with a simple 'center-of-mass' is if there are particles from both sides of the box
    x0 = np.median(x)
    y0 = np.median(y)
    z0 = np.median(z)
    cen = np.array([x0, y0 ,z0])
    
    dx = periodic_dx(x,x0,box)
    dy = periodic_dx(y,y0,box)
    dz = periodic_dx(z,z0,box)
    
    com_x = np.sum(dx*wt)
    com_y = np.sum(dy*wt)
    com_z = np.sum(dz*wt)
    cen += np.array([com_x, com_y, com_z]);

    return cen
    


## simple function to wrap x/y/z coordinates near the box edges, to give 
##   the 'correct' x-x0 in periodic coordinates
def periodic_dx( x, x0, box ):
    b0 = box / 2.
    dx = x - x0
    too_large = (dx > b0)
    dx[too_large] -= box
    too_small = (dx < -b0)
    dx[too_small] += box
    return dx
    
    

def get_galaxy_properties(PPP_type, PPP_header, ptype, cen, r_cut_halo, \
        box_size=1.0, N_in_halo_min = 10 , r_eff_cut=1.e40 ):
    
    dx = periodic_dx( PPP_type['p'][:,0], cen[0], box_size )
    dy = periodic_dx( PPP_type['p'][:,1], cen[1], box_size )
    dz = periodic_dx( PPP_type['p'][:,2], cen[2], box_size )
    r2 = dx*dx + dy*dy + dz*dz

    r2_cut = r_cut_halo*r_cut_halo
    in_halo = (r2 < r2_cut)
    N_in_halo = r2[in_halo].size

    ## set dummy/zero values for various quantities so there is something to pass
    sfr = 0.
    sf_time = 0.
    rho = 0.
    r_half = 0.
    u_internal = 0.
    vel_disp = 0.
    sf_time_rhalf = 0.
    r_median = 0.
    r_rms = 0.
    m_in_halo = 0.
    m_reff = 0.
    p_com = np.zeros(3)
    v_com = np.zeros(3)
    j_angmom = np.zeros(3)
    frac_time = np.array([0.95, 0.9, 0.7, 0.5, 0.25, 0.1])
    time_cuts = 0.*frac_time
    m_formed_since = 0.*frac_time
    
    n_species=0; 
    if (len(PPP_type['z'].shape) > 1): n_species = PPP_type['z'].shape[1]
    if (n_species<=1):
        z_metal = z_metal_rhalf = 0.0
    else:
        z_metal = z_metal_rhalf = np.zeros(n_species)

    if (N_in_halo > N_in_halo_min):
        ## redo the centering since the baryons may not be in the same place as the DM:
        cen_mod = recalculate_zoom_center( dx, dy, dz, cen=[0.,0.,0.], clip_size=r_cut_halo )
        cen += cen_mod; dx-=cen_mod[0]; dy-=cen_mod[1]; dz-=cen_mod[2];
        r2_all=dx*dx+dy*dy+dz*dz; r2=r2_all[in_halo]
        m = PPP_type['m'][in_halo]
        m_in_halo = np.sum(m)

        m_wt = m / m_in_halo
        ## sort in radius
        sort_r = np.argsort(r2)

        ## sizes / scale length(s) : 
        r_median = np.sqrt( np.median(r2) ) 
        r_rms = np.sqrt( np.sum( m_wt * r2 ) )
        m_cumulative = np.cumsum( m[sort_r] )
        f_cumulative = m_cumulative / m_in_halo
        r2_sort=r2[sort_r]
        #r_half = np.sqrt( np.interp( 0.5, f_cumulative, r2_sort ) ) ## true half-mass radius
        r_half = np.sqrt( np.min( r2_sort[f_cumulative >= 0.5 ] ) )

        ## define a weight function that doesn't evaluate outside the 1/2 scale-length, 
        ##   useful for looking at gradients and/or whether there is a difference in small-r
        in_halo_rhalf = (r2_all < r_half*r_half)
        m_wt_rhalf = PPP_type['m'][in_halo_rhalf]
        m_wt_rhalf /= np.sum(m_wt_rhalf)
        in_reff_cut = (r2_all < r_eff_cut*r_eff_cut) & (r2_all < r2_cut)

        ## metallicity
        if ((ptype==0) | (ptype==4)):
            if (len(PPP_type['z'].shape) > 1):
                ## should correctly sum for multiple species
                n_species = PPP_type['z'].shape[1]
                q = PPP_type['z'][in_halo,:]
                q_rhalf = PPP_type['z'][in_halo_rhalf,:]
                z_metal = np.zeros(n_species)
                z_metal_rhalf = np.zeros(n_species)
                for j_z in range(n_species):
                    z_metal[j_z] = np.sum( m_wt * q[:,j_z], axis=0 )
                    z_metal_rhalf[j_z] = np.sum( m_wt_rhalf * q_rhalf[:,j_z], axis=0 )
            else:
                z_metal = np.sum( m_wt * PPP_type['z'][in_halo], axis=0 )
                z_meta_rhalf = np.sum( m_wt_rhalf * PPP_type['z'][in_halo_rhalf], axis=0 )
    
        if (ptype==0):
            ## star formation rate
            sfr = np.sum( PPP_type['sfr'][in_halo] )
            ## density, temperature, etc ? 
        
        if (ptype==4):
            ## stellar formation times
            sf_time = np.median( PPP_type['age'][in_halo] )
            sf_time_rhalf = np.sum( m_wt_rhalf * PPP_type['age'][in_halo_rhalf] )
            # sf_recent = 
            # (some cut on stellar mass formed in some short delta_t):
            time = PPP_header['time']
            age = PPP_type['age'][in_halo]
            time_cuts = time*frac_time
            for j,frac_time_ok in zip(range(frac_time.size),frac_time):
                ok = (age < time) & (age > time*frac_time_ok)
                m_formed_since[j] = np.sum(m[ok])
    
        ## kinematics 
        ##   (here, its useful to cut to inside/around r_half)
        m = PPP_type['m'][in_halo_rhalf]
        m_reff = np.sum( PPP_type['m'][in_reff_cut] )

        ## center of mass of baryons (be careful of box-wrapping!)
        pos_0 = PPP_type['p'][in_halo_rhalf]
        dx = periodic_dx( pos_0[:,0], cen[0], box_size )
        dy = periodic_dx( pos_0[:,1], cen[1], box_size )
        dz = periodic_dx( pos_0[:,2], cen[2], box_size )
        p_com = np.zeros(3)
        p_com[0] = np.sum( m_wt_rhalf * dx )
        p_com[1] = np.sum( m_wt_rhalf * dy )
        p_com[2] = np.sum( m_wt_rhalf * dz )
        pos = np.zeros([dx.size,3]); 
        pos[:,0]=dx-p_com[0]; pos[:,1]=dy-p_com[1]; pos[:,2]=dz-p_com[2];
        p_com += cen       

        N_in_half = m.size        
        v_com=np.zeros(3)
        vel = PPP_type['v'][in_halo_rhalf]
        ## velocity centroid and dispersion (much simpler since box-wrapping doesn't matter)
        for j in [0,1,2]: 
            v_com[j] = np.sum( m_wt_rhalf * vel[:,j] )
            vel[:,j] -= v_com[j]
        v2 = vel[:,0]*vel[:,0] + vel[:,1]*vel[:,1] + vel[:,2]*vel[:,2]
    
        ## velocity dispersion
        vel_disp = np.sqrt( np.sum( m_wt_rhalf * v2 ) )
    
        ## specific angular momentum 
        j_spec = np.cross( pos, vel, axis=1 )
        for j in [0,1,2]:
            j_angmom[j] = np.sum( m_wt_rhalf * j_spec[:,j] )
        
    ## return everything as a dictionary:
    return {\
        'n':N_in_halo,'m':m_in_halo,'m_reff':m_reff,'r_half':r_half,'v_disp':vel_disp,\
        'pos':p_com,'vel':v_com,'j':j_angmom,'z':z_metal,'z_half':z_metal_rhalf,\
        'sfr':sfr,'age':sf_time,'age_half':sf_time_rhalf,'time_cuts':time_cuts,\
        'mass_per_time_cut':m_formed_since}



def recalculate_zoom_center( x0, y0, z0, cen=[0.,0.,0.], clip_size=2.e10 ):
    #rgrid=np.array([1.0e10,1000.,700.,500.,300.,200.,100.,70.,50.,30.,20.,10.,5.,2.5,1.]); jmax_loop=5;
    rgrid=np.array([1.0e10,1000.,500.,250.,100.,50.,25.,10.,5.,2.5,1.]); jmax_loop=2;
    rgrid = rgrid[rgrid <= clip_size];
    r2grid = rgrid * rgrid;
    for i_rcut in range(r2grid.size):
        for j_looper in range(jmax_loop):
            x=x0-cen[0]; y=y0-cen[1]; z=z0-cen[2];
            r2 = x*x + y*y + z*z
            ok = (r2 < r2grid[i_rcut])
            if (r2[ok].size > 1000):
                x=x[ok]; y=y[ok]; z=z[ok];
                if (i_rcut <= (rgrid.size-5)):
                    cen += np.array([np.median(x),np.median(y),np.median(z)]);
                else:
                    cen += np.array([np.mean(x),np.mean(y),np.mean(z)]);
            else:
                if (r2[ok].size > 10):
                    x=x[ok]; y=y[ok]; z=z[ok];
                    cen += np.array([np.mean(x),np.mean(y),np.mean(z)]);
    return cen;



def load_halo_group_ids(snapdir,snapnum, \
        cosmological=1,h0=1,skip_bh=1,four_char=0):

    ## first check if we already have this saved:
    fname = get_halo_id_filename(snapdir, snapnum, four_char=four_char)
    if os.path.exists(fname):
        ## file exists, load it up
        halo_id_file = h5py.File(fname,'r')
        halo_ids = np.copy(np.array(halo_id_file["HaloIDs"]))
        halo_id_file.close()
    else:
        ## make them ourselves
        halo_ids = dump_halo_finder_results(snapdir,snapnum, \
            cosmological=cosmological,h0=h0,skip_bh=skip_bh,four_char=four_char)
    
    return halo_ids



def dump_halo_finder_results(snapdir,snapnum, \
        cosmological=1,h0=1,skip_bh=1,four_char=0):
    
    ## load the relevant snapshot: collating the DM particles:
    have_key = 0
    for ptype in [1,2,3,5]:
        P=gadget.readsnap(snapdir,snapnum,ptype,\
            h0=h0,cosmological=cosmological,skip_bh=skip_bh,four_char=four_char)
        if (P['k']==1): 
            if (ptype<5) or (P['m'].shape[0]>1000):
                if (have_key==0):
                    pos=P['p']; m=P['m']; have_key=1;
                else:
                    pos=np.concatenate((pos,P['p'])); 
                    m=np.concatenate((m,P['m']));

    if (pos.size <= 1): 
        print('Not enough particles! only n=',pos.size)
        return -1
    
    ## renormalize quantities (for halo finder convenience)
    pmin = np.min(pos);
    pmax = np.max(pos);
    pos = (pos-pmin)/(pmax-pmin);
    m /= np.sum(m);
    x=pos[:,0]; y=pos[:,1]; z=pos[:,2];

    ## now call the actual halo finder and get group membership info for DM particles:
    particle_parent_group_ids = call_halo_finder(x,y,z,m, threshold_for_halo=160.)

    ## and write this out to a saved file
    fname = get_halo_id_filename(snapdir, snapnum, four_char=four_char)
    outfi = h5py.File(fname,'w')
    dset_im = outfi.create_dataset('HaloIDs',data=particle_parent_group_ids)
    outfi.close()

    return particle_parent_group_ids
    


def call_halo_finder( x, y, z, mass, threshold_for_halo = 160.0):
    ## load the routine we need
    exec_call=util.return_python_routines_cdir()+'/hop_halofinder/hop.so'
    routine=ctypes.cdll[exec_call];
    ## ctypes is very picky about its data types, 
    ##    so treat the following to pre-process appropriately:
    x=fcor(x); y=fcor(y); z=fcor(z); m=fcor(mass);
    N=np.int(x.size); threshold_for_halo=np.float(threshold_for_halo);
    print(N)
    print(len(mass))
    ## cast the variables to store the results
    out_cast = ctypes.c_int*N; 
    particle_parent_group_ids_out = out_cast();
    ## main call to the calculation routine
    routine.hop_shared( ctypes.c_int(N), \
        vfloat(x), vfloat(y), vfloat(z), vfloat(m), \
        ctypes.byref(particle_parent_group_ids_out), \
        ctypes.c_float(threshold_for_halo) );
    ## now put the output arrays into a useful format 
    particle_parent_group_ids = np.copy(np.ctypeslib.as_array(particle_parent_group_ids_out));
    return particle_parent_group_ids


def test_fof():
    sdir='/Users/phopkins/Downloads/m12_zoom/z5_box'
    snum=29
    ptypes=[4]
    linklength=1.e-5 # ~1pc for snap29 of z5_box
    linklength=3.e-5
    #ids=load_fof_group_ids(sdir,snum,ptypes,linklength,cosmological=1)
    compile_fof_group_properties(sdir,snum,ptypes,linklength,cosmological=1)


def compile_fof_group_properties(snapdir, snapnum, ptypes, linking_length, \
        cosmological=1,h0=1,skip_bh=1,four_char=0):

    ## first load the fof ID set (or build it, if none exists): 
    halo_ids_dm,linking_length = load_fof_group_ids(snapdir,snapnum,ptypes,linking_length, \
            cosmological=cosmological,h0=h0,skip_bh=skip_bh,four_char=four_char)
    
    ## load the for the snapshot (so we don't have to re-do this in-loop):
    ##   -- collate particles of the type we care about -- 
    have_key = 0
    j_check=0;
    while (have_key==0):
        P_header = gadget.readsnap(snapdir,snapnum,j_check,header_only=1,\
            h0=h0,cosmological=cosmological,skip_bh=skip_bh,four_char=four_char)
        if (P_header['k'] != -1): have_key=1;
        j_check += 1;
        if (j_check > 5): have_key=1;
    have_key = 0
    for ptype in ptypes:
        P=gadget.readsnap(snapdir,snapnum,ptype,\
            h0=h0,cosmological=cosmological,skip_bh=skip_bh,four_char=four_char)
        if (P['k']==1): 
            if (have_key==0):
                pos=P['p']; m=P['m']; type=np.zeros(P['m'].size)+np.float(ptype)
                have_key = 1
            else:
                pos=np.concatenate((pos,P['p'])); 
                m=np.concatenate((m,P['m']));
                type=np.concatenate((type,np.zeros(P['m'].size)+np.float(ptype)));
    m_dm=m; x_dm=pos[:,0]; y_dm=pos[:,1]; z_dm=pos[:,2]; type_dm=type;
    if(cosmological==1):
        box_size=P_header["boxsize"]
    else:
        box_size=1.e20;
    print('box_size = ',box_size)

    ## ok now loop through and flag everything we care about in the files:
    min_number_of_particles_in_halo = 30

    ## prepare the dictionaries::
    hhh_dm = np.empty(0,dtype=object)
    halo_id_max = max(halo_ids_dm)
    time = P_header['time']
    print('N_GROUPS == ',halo_id_max)
    print('.. beginning loop now...')
    for halo_id in range(halo_id_max+1):
        in_halo = ( halo_ids_dm == halo_id )
        m = m_dm[in_halo]
    
        ## halos are already sorted by total particle number, so once we get to those 
        ##   which are 'too small', we can terminate this loop
        N_dm_halo = m.size
        if (N_dm_halo < min_number_of_particles_in_halo): continue
        m_dm_halo = np.sum(m)
        print(halo_id, min_number_of_particles_in_halo, m.size)
        
        ## now get the halo position
        x = x_dm[in_halo]; y = y_dm[in_halo]; z = z_dm[in_halo]; wt = m/m_dm_halo;
        cen = get_halo_com(x, y, z, wt, box_size)
    
        dx=periodic_dx(x,cen[0],box_size)
        dy=periodic_dx(y,cen[1],box_size)
        dz=periodic_dx(z,cen[2],box_size)
        r2 = dx*dx + dy*dy + dz*dz
        rmax = np.sqrt(np.max(r2))
        rrms = np.sqrt(np.sum(m*r2)/m_dm_halo)
        rmed = np.sqrt(np.median(r2))
        ## not the same as half-mass, if multi-mass particles
    
        ## now loop over the groups and evaluate the relevant properties
        
        ggg_dm = {'id':halo_id,'n':N_dm_halo,'m':m_dm_halo,'pos':cen,'r_max':rmax,\
            'r_med':rmed,'r_rms':rrms,'time':time}
        print('(( ',halo_id,' ))')
        print(ggg_dm)
        hhh_dm = np.append( hhh_dm, ggg_dm )
        
        print(' ')
        print(' ')
        sys.stdout.flush()

    ## ok done with major loop over halos
    N_halos_in_final = hhh_dm.size
    print('... done with census: N_groups_final = ',N_halos_in_final)

    if (2==0):
        ## open the file to record all saved properties (so they can be easily loaded later)
        fname = get_halo_id_filename(snapdir,snapnum,full_halo_results=1,four_char=four_char)
        print('writing to file : ',fname)
        outfi = h5py.File(fname,'w')
    
        ## organize keywords :: 
        type_vec=[0,1,2]
        if (dm_only==1): type_vec=[0]
        for type_loop in type_vec:
            if (type_loop==0):
                dic = hhh_dm
                prefix = 'halo_'
            if (type_loop==1):
                dic = hhh_gas
                prefix = 'gas_'
            if (type_loop==2):
                dic = hhh_star
                prefix = 'star_'
            for key_var in dic[0].keys():
                ## fill in an array with this keyname
                d00 = np.empty(0)
                for i in range(N_halos_in_final): 
                    d0_tmp = dic[i]
                    d00 = np.append(d00,d0_tmp[key_var])
                d0 = np.array(d0_tmp[key_var])
                if (d0.size>1): d00=d00.reshape(-1,d0.size)
                dataset = outfi.create_dataset(prefix+key_var,data=d00)
        outfi.close() ## done writing
    print('done!')

    ## (to read) :::
    #infiname = get_halo_id_filename(snapdir,snapnum,full_halo_results=1,four_char=four_char)
    #infi=h5py.File(infiname,'r')
    ## now any keyword can be called by its name
    #sigma = np.array(infi["star_v_disp"])   # etc
    ## use 'list(infi)' to see available keys
    #infi.close()
    
    return 1



def load_fof_group_ids(snapdir,snapnum,ptypes,linking_length, \
        cosmological=1,h0=1,skip_bh=1,four_char=0):

    ## first check if we already have this saved:
    ext='_fof_p'
    for ptype in ptypes: ext=ext+str(ptype);
    fname = get_halo_id_filename(snapdir, snapnum, four_char=four_char) + ext
    if os.path.exists(fname):
        ## file exists, load it up
        halo_id_file = h5py.File(fname,'r')
        fof_ids = np.copy(np.array(halo_id_file["fof_ids"]))
        linking_length = np.copy(np.array(halo_id_file["fof_linkinglength"]))
        halo_id_file.close()
    else:
        ## make them ourselves
        fof_ids = dump_fof_finder_results(snapdir,snapnum,ptypes,\
            linking_length=linking_length, \
            cosmological=cosmological,h0=h0,skip_bh=skip_bh,four_char=four_char)
    
    return fof_ids, linking_length


def dump_fof_finder_results(snapdir,snapnum,ptypes, \
        linking_length=0.2, \
        cosmological=1,h0=1,skip_bh=1,four_char=0):
    
    ## load the relevant snapshot: collating the particles of interest:
    have_key = 0
    ext='_fof_p'
    for ptype in ptypes:
        P=gadget.readsnap(snapdir,snapnum,ptype,\
            h0=h0,cosmological=cosmological,skip_bh=skip_bh,four_char=four_char)
        if (P['k']==1): 
            ext = ext+str(ptype)
            if (have_key==0):
                pos=P['p']; m=P['m']; have_key=1;
            else:
                pos=np.concatenate((pos,P['p'])); 
                m=np.concatenate((m,P['m']));

    if ((have_key==0)|(pos.size <= 1)): 
        print('Not enough particles! only n=',pos.size)
        return -1
    
    ## renormalize quantities (for group finder convenience)
    pmin = np.min(pos);
    pmax = np.max(pos);
    print('pMAXMIN == ',pmax,pmin,pmax-pmin)
    pos = (pos-pmin)/(pmax-pmin);
    m /= np.sum(m);
    x=pos[:,0]; y=pos[:,1]; z=pos[:,2];

    ## now call the actual halo finder and get group membership info:
    particle_parent_group_ids = call_fof_finder(x,y,z,m,linking_length=linking_length)

    ## and write this out to a saved file
    fname = get_halo_id_filename(snapdir, snapnum, four_char=four_char) + ext
    outfi = h5py.File(fname,'w')
    dset_im1 = outfi.create_dataset('fof_ptype',data=ptypes)
    dset_im2 = outfi.create_dataset('fof_linkinglength',data=linking_length)
    dset_im3 = outfi.create_dataset('fof_ids',data=particle_parent_group_ids)
    outfi.close()

    return particle_parent_group_ids
    
    


def call_fof_finder( x, y, z, mass, linking_length = 0.2):
    ## load the routine we need
    exec_call=util.return_python_routines_cdir()+'/fof/fof.so'
    routine=ctypes.cdll[exec_call];
    ## ctypes is very picky about its data types, 
    ##    so treat the following to pre-process appropriately:
    x=fcor(x); y=fcor(y); z=fcor(z); m=fcor(mass);
    N=np.int(x.size); linking_length=np.float(linking_length);
    ## cast the variables to store the results
    out_cast = ctypes.c_int*N; 
    particle_parent_group_ids_out = out_cast();
    ## main call to the calculation routine
    routine.fof_shared( ctypes.c_int(N), \
        vfloat(x), vfloat(y), vfloat(z), vfloat(m), \
        ctypes.byref(particle_parent_group_ids_out), \
        ctypes.c_float(linking_length) );
    ## now put the output arrays into a useful format 
    particle_parent_group_ids = np.copy(np.ctypeslib.as_array(particle_parent_group_ids_out));
    return particle_parent_group_ids



def fcor(x):
    return np.array(x,dtype='f',ndmin=1)

def vfloat(x):
    return x.ctypes.data_as(ctypes.POINTER(ctypes.c_float));

def snap_ext(snum,four_char=0):
	ext='00'+str(snum);
	if (snum>=10): ext='0'+str(snum)
	if (snum>=100): ext=str(snum)
	if (four_char==1): ext='0'+ext
	if (snum>=1000): ext=str(snum)
	return ext;

def get_halo_id_filename(snapdir, snum, full_halo_results=0, four_char=0):
    ## (first parse the directory names (for file-naming conventions))
    s0=snapdir.split("/"); snapdir_specific=s0[len(s0)-1]; n_s0=1;
    if(len(snapdir_specific)<=1): snapdir_specific=s0[len(s0)-2]; n_s0=2;
    ss=snap_ext(snum,four_char=four_char)
    fname_base=snapdir+'/'+snapdir_specific+ss
    
    ext = '.dmhalo_ids'
    if (full_halo_results==1): ext = '.dmhalo_prop'
    
    return fname_base+ext

