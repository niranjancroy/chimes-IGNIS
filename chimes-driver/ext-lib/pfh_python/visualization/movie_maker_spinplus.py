import numpy as np
import matplotlib
import ctypes
import math
import h5py
import os.path
import struct
import array
import sys
import gc
from subprocess import call
import matplotlib.colors
import matplotlib.cm
import matplotlib.pyplot as plt
import visualization.colors as viscolors
import visualization.get_attenuated_stellar_luminosities as getlum
import visualization.make_threeband_image as makethreepic
import visualization.contour_makepic as cmakepic
import visualization.raytrace_projection as rayproj
import gadget_lib.load_stellar_hsml as starhsml
import pfh_utils as util
import gadget




def construct_spin_movie_single_sim(
    smaster = '/scratch/01799/phopkins/m12i_hybrid_test/m12_awetzel/m12i/fb-sym/', \
    omaster = '/work/01799/phopkins/images/', \
    sdir = 'm12i_ref13', \
    subfolder_ext = 'output', \
    snum = 600, \
    key = 'star', \
    nf_base = 100, \
    add_key = '', \
    #x00_in = 2.0, \
    x00_in = 3.0, \
    x00_med = 30.0, \
    frac_to_rotate_initial_spin = 0.76, \
    pixels = 800., \
    theta_initial = 11., \
    theta_median = 109., \
    #theta_median = 104., \
    phi_initial = 0., \
    phi_median = -60., \
    #phi_median = -55., \
    frame_min = 0, \
    frame_max = -1,\
    scattered_fraction = 0.05,\
    h_rescale_factor=1.,\
    h_max=0,
    a_scale_initial=0.90,
    a_scale_final=1.0):
    
    cen=get_snapshot_precomputed_z0_cen(sdir); print('cen = ',cen)
    if(cen[0]==0.):
        cen,n00 = fast_zoom_center(smaster+'/'+sdir+'/'+subfolder_ext+'/',snum,cosmological=1);
        print('calculated z=0 center for movie at: ',cen)

    id_i,m_i,x_i,y_i,z_i,vx_i,vy_i,vz_i,h_i,c_i,zm_i = \
        load_pos_vel_etc(smaster+'/'+sdir+'/'+subfolder_ext+'/', snum, filename_prefix = sdir, \
            h0=1,four_char=0,cosmological=1,GAS=0,BOX_MAX=15.,CENTER_BOX=cen)
    x=x_i-cen[0]; y=y_i-cen[1]; z=z_i-cen[2]; r=np.sqrt(x*x+y*y+z*z)
    ok=np.where(r < 15.)
    x=x[ok]; y=y[ok]; z=z[ok]; r=r[ok]; m=m_i[ok]; 
    vx=vx_i[ok]; vy=vy_i[ok]; vz=vz_i[ok];
    vx-=np.median(vx); vy-=np.median(vy); vz-=np.median(vz);
    jx=vy*z-vz*y; jy=vz*x-vx*z; jz=vx*y-vy*x;
    jxm=np.sum(m*jx); jym=np.sum(m*jy); jzm=np.sum(m*jz);
    jmm=np.sqrt(jxm*jxm+jym*jym+jzm*jzm)
    jxm/=jmm; jym/=jmm; jzm/=jmm;
    z_vec = np.array([jxm,jym,jzm])
    x_vec,y_vec = return_perp_vectors(z_vec)
    
    # now set up the spin
    nf1 = nf_base
    nf2 = nf_base
    nf2x = np.int(np.round(nf_base/1.5))
    nf3 = nf2
    frac1 = frac_to_rotate_initial_spin
    t0=theta_initial; t1=theta_median; p0=phi_initial; p1=phi_median;
    nf2 = np.int(np.round(nf_base*1.5))
    nf2x = np.int(np.round(nf_base*1.0))
    nf3 = nf_base
    
    n = 2.*nf1
    dt = frac1*(t1-t0)/n; dt0=dt;
    dp = frac1*(p1-p0)/n; dp0=dp;
    dlnx = 2.5/n; dlnx0=dlnx;
    theta_grid = t0 + dt*np.arange(0,n)
    phi_grid = p0 + dp*np.arange(0,n)
    x00_grid = x00_med * np.exp(dlnx*(n-np.arange(0,n)))

    n = nf2
    dt = 1.*(1.-frac1)*(t1-t0)/n;
    dp = (1.-frac1)*(p1-p0)/n;
    theta_grid,dt0 = tcompile(theta_grid,dt,dt0,n)
    phi_grid,dp0 = tcompile(phi_grid,dp,dp0,n)
    x00_grid,dlnx0 = xcompile(x00_grid,x00_in,dlnx0,n)

    n = nf2x
    dt = 1.*(1.-frac1)*(t1-t0)/n
    dp = 2.*(1.-frac1)*(p1-p0)/n
    theta_grid,dt0 = tcompile(theta_grid,dt,dt0,n)
    phi_grid,dp0 = tcompile(phi_grid,dp,dp0,n)
    x00_grid,dlnx0 = xcompile(x00_grid,0.5*x00_med,dlnx0,n)

    n = nf1
    dt = (frac1-1.*(1.-frac1))*frac1*(t1-t0)/n
    dp = (frac1-2.*(1.-frac1))*(p1-p0)/n
    theta_grid,dt0 = tcompile(theta_grid,dt,dt0,n)
    phi_grid,dp0 = tcompile(phi_grid,dp,dp0,n)
    x00_grid,dlnx0 = xcompile(x00_grid,0.8*x00_med,dlnx0,n)

    n = nf3
    dt = 0.33*frac1*(t1-t0)/n
    dp = 0.5*180.0/n
    xf = 1.2

    n *= 2.5
    dt *= 1.5
    dp *= 1.5
    xf *= 5.

    theta_grid,dt0 = tcompile(theta_grid,dt,dt0,n)
    phi_grid,dp0 = tcompile(phi_grid,dp,dp0,n)
    x00_grid,dlnx0 = xcompile(x00_grid,xf*x00_med,dlnx0,n)
    
    qq = 1.*np.arange(0,x00_grid.size)
    a_grid = a_scale_initial + (a_scale_final-a_scale_initial) * qq/np.max(qq)

    print('ESTIMATED_NUMBER_OF_FRAMES == ',x00_grid.size)

    addlayer='';set_added_layer_alpha=0.0;lightk=0;dustgas=1;ctable='heat_purple'
    threecolor=1; invert_colors=0;
    if(key=='star'):
        dynr = 2.e3 + 0.*x00_grid; 
        maxd = 1.1e-2 * (x00_grid / 30.)**(-1.5); 
        dustgas=2.2; lightk=0; 
        
        if (add_key=='nodust'): 
            dustgas = 1.e-5
        if (add_key=='CO'): 
            addlayer='CO'; ctable='heat_blue'; set_added_layer_alpha=0.4;
        if (add_key=='Halpha'): 
            addlayer='Halpha'; ctable='heat_green'; set_added_layer_alpha=0.3;
        if (add_key=='SFR'):
            addlayer='SFR'; ctable='heat_purple'; set_added_layer_alpha=0.3;
        if (add_key=='Xray'):
            addlayer='Xray'; ctable='heat_redyellow'; set_added_layer_alpha=0.4;
        if (add_key=='Zmetal'):
            addlayer='Zmetal'; ctable='heat_orange'; set_added_layer_alpha=0.8;
            
    if(key=='gas'):
        dynr = 1.e3 + 0.*x00_grid; 
        maxd = 6.0e-2 * (x00_grid / 30.)**(-1.5);
        threecolor=0; maxd*=10.0; dynr*=3.; lightk=0; invert_colors=0;

    movie_maker(frame_min=frame_min,frame_max=frame_max,
        snapdir=smaster+'/'+sdir,
        outputdir_master=omaster,pixels=pixels,subfolder_ext=subfolder_ext,
        do_with_colors=1,threecolor=threecolor,cosmological=1,
        show_gasstarxray=key,include_lighting=lightk,
        set_added_layer_ctable=ctable,add_gas_layer_to_image=addlayer,
        dust_to_gas_ratio_rescale=dustgas,set_added_layer_alpha=set_added_layer_alpha,
        invert_colors=invert_colors,
        show_scale_label=0,show_time_label=0,spin_movie=1,
        scattered_fraction=scattered_fraction, 
        h_rescale_factor=h_rescale_factor,h_max=h_max,
        projection_vector_x=x_vec,projection_vector_y=y_vec,projection_vector_z=z_vec,
        a_grid_topass=a_grid,dx_box_grid_topass=x00_grid,
        theta_grid_topass=theta_grid,phi_grid_topass=phi_grid,
        maxden_grid_topass=maxd,dynrange_grid_topass=dynr)


##
## here is the main routine for the galaxy simulation movies: 
##   it builds the snapshot list, does the centering, goes through frames and calls
##   the appropriate routines for interpolation between snapshots, and then calls the 
##   image_maker routine to render the frame: many switches here (most don't need to be set), 
##   but tune it as you need to (I recommend a low snapshot spacing to test)
##


def movie_maker(\
    xmax_of_box=50., #"""scale (in code units) of movie box size"""
    snapdir='/scratch/01799/phopkins/m12i_hybrid_test/m12_awetzel/m12i/fb-sym/m12i_ref12', #"""location of snapshots"""
    outputdir_master='/work/01799/phopkins/images/', #"""parent folder for frames dump"""
    subfolder_ext='output', #"""subfolder of snapdir containing the snapshots"""
    show_gasstarxray='stars', #"""determines if image is 'stars','gas','xr' etc"""
    frames_per_gyr=100., #"""sets spacing of frames (even in time)"""
    cosmological=1, #"""is this is a cosmological (comoving) simulation?"""
    time_min=0, #"""initial time of movie"""
    time_max=0, #"""final time (if<=0, just set to time of maximum snapshot)"""
    frame_min=0, #"""initial frame number to process"""
    frame_max=1.0e10, #"""final frame number to process"""
    i_snap_min=0, #"""minimum snapshot number in list to scan"""
    i_snap_max=0, #"""maximum snapshot number to scan (if =0, will just use largest in snapdir)"""
    pixels=720, #"""images will be pixels*pixels"""
    show_time_label=0, #"""do or don't place time label on each frame"""
    show_scale_label=0, #"""do or don't place physical scale label on each frame"""
    temp_max=1.0e6, #"""maximum gas temperature for color-weighted gas images (as 'temp_cuts' in 3-color)"""
    temp_min=3.0e2, #"""minimum gas temperature for color-weighted gas images (as 'temp_cuts' in 3-color)"""
    theta_0=90., #"""(initial) polar angle (if rotating)"""
    phi_0=90., #"""(initial) azimuthal angle (if rotating)"""
    center_on_bh=0, #"""if have a bh particle, uses it for the centering"""
    center_on_com=0, #"""simple centering on center of mass"""
    set_fixed_center=0, #"""option to set a fixed center for the whole movie"""
    sdss_colors=0, #"""use sdss color-scheme for three-color images"""
    nasa_colors=1, #"""use nasa color-scheme for three-color images"""
    use_h0=1, #"""correct snapshots to physical units"""
    four_char=0, #"""snapshots have 4-character filenames"""
    skip_bh=0, #"""skips reading bh information in snapshots"""
    use_old_extinction_routine=0, #"""uses older (faster but less stable) stellar extinction routine"""
    add_extension='', #"""extension for movie snapshot images"""
    do_with_colors=1, #"""make color images"""
    threecolor=1, #"""do or dont use the standard three-color projection for images"""
    camera_opening_angle=45., #"""camera opening angle (specifies camera distance)"""
    scattered_fraction=0.01, #"""fraction of light that is unattenuated"""
    z_to_add=0.0, #"""add this metallicity to all particles"""
    #min_stellar_age=0.0, #"""force star particles to be at least this old to avoid interpolation artifacts"""
    min_stellar_age=1.e-3, #"""force star particles to be at least this old to avoid interpolation artifacts"""
    include_lighting=0, dust_to_gas_ratio_rescale=1.,invert_colors=0,spin_movie=0,h_rescale_factor=1.,h_max=0, #"""stuff"""
    add_gas_layer_to_image='', set_added_layer_alpha=0.3, set_added_layer_ctable='heat_purple', #"""stuff"""
    projection_vector_x=[1,0,0],projection_vector_y=[0,1,0],projection_vector_z=[0,0,1], #"""stuff"""
    a_grid_topass=[0],t_grid_topass=[0],dx_box_grid_topass=[0],theta_grid_topass=[0],phi_grid_topass=[0],maxden_grid_topass=[0],dynrange_grid_topass=[0], #"""stuff"""
    use_polar_interpolation=1):

    ## first decide whether its a gas or stellar image
    SHOW_STARS=0; SHOW_GAS=1;
    if((show_gasstarxray=='star') or (show_gasstarxray=='stars') or (show_gasstarxray=='st')): 
        SHOW_STARS=1; SHOW_GAS=0;

    ## parse the directory names (for file-naming conventions)
    s0=snapdir.split("/"); 
    snapdir_specific=s0[len(s0)-1]; n_s0=1; 
    if(len(snapdir_specific)<=1): snapdir_specific=s0[len(s0)-2]; n_s0=2;
    snapdir_master=''; 
    for j in s0[0:len(s0)+1-n_s0]: snapdir_master += str(j)+'/';
    outputdir_master+='/' ## just to be safe
    full_snapdir = snapdir+'/'+subfolder_ext
    
    ## build the snapshot list in the directory of interest
    print('... building snapshot list in directory ...')
    snapshot_list = build_snapshot_list(full_snapdir,four_char=four_char )
    if(i_snap_max<=0): i_snap_max=snapshot_list[snapshot_list.size-2];

    ## to avoid 'camera jitter', the safest thing (albeit more expensive) is to 
    ##   pre-compute the movie center at all snapshots, then use a smoothed tracking of it
    ##   -- call subroutine to check if its done for us, and if not, build it --
    print('... building/loading camera positions list ...')
    build_camera_centering( snapshot_list,full_snapdir,snapdir_specific,outputdir_master,
        force_rebuild=0,four_char=four_char,cosmological=cosmological)

    a_grid_topass = np.array(a_grid_topass)
    t_grid_topass = np.array(t_grid_topass)
    if(a_grid_topass.size+t_grid_topass.size > 2):
        #ok = (a_grid_topass >= time_min)&(a_grid_topass <= time_max)
        ok = (a_grid_topass >= -1)
        a_scale_grid = a_grid_topass[ok]
        time_frame_grid = cosmological_time(a_scale_grid)
        t = 0.0*time_frame_grid
        x_plot_scale_grid = np.array(dx_box_grid_topass)[ok]
        y_plot_scale_grid = 1.*x_plot_scale_grid
        theta_frame_grid = np.array(theta_grid_topass)[ok]
        phi_frame_grid = np.array(phi_grid_topass)[ok]
        maxden_grid = np.array(maxden_grid_topass)[ok]
        dynrange_grid = np.array(dynrange_grid_topass)[ok]
    else:
        print('... building list of times for movie frames ...')
        ## now build the times for each frame of the movie
        time_frame_grid, a_scale_grid = build_time_frame_grid( snapdir, snapshot_list, \
              frames_per_gyr, time_min=time_min, time_max=time_max, \
              use_h0=use_h0,four_char=four_char,cosmological=cosmological,skip_bh=skip_bh)

        ##
        ## some (simple) defaults for the scaling of dynamic range, etc
        ##  (can do fancy things for movie panning/zooming etc by changing the variables here)
        ##
        t=0.0*time_frame_grid
        x_plot_scale_grid = t + xmax_of_box
        y_plot_scale_grid = x_plot_scale_grid
        #theta_frame_grid = t + 90. # should be in degrees
        theta_frame_grid = t + theta_0 # should be in degrees
        #phi_frame_grid = t + 90. ## x-y
        phi_frame_grid = t + phi_0 ## x-y
        dynrange_0 = 3.0e2
        maxden_0 = 3.8e7
        maxden_grid = t + maxden_0 * ((50./x_plot_scale_grid)**(0.3))
        if(SHOW_STARS==1):
            dynrange_0 *= 100.0 #XXX CCH: originally 3.0
            maxden_grid *= 100.0 #XXX CCH: originally 3.0
        maxden_grid /= 1.0e10 # recalling in gadget units, m=1.0e10 m_sun
        # layer below is for zoom_dw large-scale GAS sims ONLY
        #z_grid = 1./a_scale_grid - 1.
        #dynrange_0 *= 20. * (1.+np.exp(-(z_grid-2.)))
        #maxden_grid *= 3.
        dynrange_grid=dynrange_0 + 0.*x_plot_scale_grid
    
    
    print('... entering main snapshot loop ...')
    ## now enter main loop over snapshots in the simulation set
    snapshot_f_number_that_i_just_loaded = -1;
    snaps_to_do = (snapshot_list >= i_snap_min) & (snapshot_list <= i_snap_max)
    ii = np.arange(snapshot_list.size)
    snap_grid_to_loop = ii[snaps_to_do]
    for i_snap in snap_grid_to_loop:
        ## load header info for snapshot 'i'
        PPP_head = gadget.readsnap(snapdir+'/'+subfolder_ext,snapshot_list[i_snap],1,header_only=1,\
            h0=use_h0,four_char=four_char,cosmological=cosmological,skip_bh=skip_bh)
        if(cosmological==1):
            a_scale_i=PPP_head['time']; t_i=cosmological_time(a_scale_i); 
        else: 
            t_i=PPP_head['time']; a_scale_i=t_i
        a_scale_i=np.array(a_scale_i); t_i=np.array(t_i);
            
        ## load header info for snapshot 'f=i+1'
        PPP_head = gadget.readsnap(snapdir+'/'+subfolder_ext,snapshot_list[i_snap+1],1,header_only=1,\
            h0=use_h0,four_char=four_char,cosmological=cosmological,skip_bh=skip_bh)
        if(cosmological==1):
            a_scale_f=PPP_head['time']; t_f=cosmological_time(a_scale_f); 
        else: 
            t_f=PPP_head['time']; a_scale_f=t_f
        a_scale_f=np.array(a_scale_f); t_f=np.array(t_f);
        
        print('... snapshot headers loaded: getting timestep info: ...')
        ## now calculate whether there are any frames in the time-range between the snapshots
        delta_t = t_f - t_i
        print('Timestep is ti/tf/delta_t ',t_i,t_f,delta_t)
        i_frame = np.array(range(time_frame_grid.size));
        check_frames_to_do = (time_frame_grid >= t_i) & (time_frame_grid <= t_f) & \
            (i_frame <= frame_max) & (i_frame >= frame_min);
        frames_to_do = i_frame[check_frames_to_do];
        n_frames_to_do = frames_to_do.size;

        if(n_frames_to_do>0):
            print('... processing ',n_frames_to_do,' frames in this snapshot interval ...')

            center_i = get_precalc_zoom_center(snapdir+'/'+subfolder_ext, snapdir_specific, outputdir_master,a_scale_i,cosmological=cosmological)
            center_f = get_precalc_zoom_center(snapdir+'/'+subfolder_ext, snapdir_specific, outputdir_master,a_scale_f,cosmological=cosmological)
            ## correct to comoving coordinates for interpolation (so correctly capture hubble flow)
            if (cosmological==1):
                center_i /= a_scale_i
                center_f /= a_scale_f
            print('... ... centered (i/f) at ',center_i,center_f)

            print('... ... loading snapshot bricks ... ...')
            ## tolerance for keeping particles outside the box, for this calculation
            #xmax_tmp=xmax_of_box*1.5 + 10.
            xmax_tmp=np.max(x_plot_scale_grid[frames_to_do])*1.5 + 10.
            ## load the actual --data-- for snapshot (i): 
            id_i,m_i,x_i,y_i,z_i,vx_i,vy_i,vz_i,h_i,c_i,zm_i = \
                load_pos_vel_etc(snapdir+'/'+subfolder_ext, snapshot_list[i_snap], filename_prefix = snapdir_specific, \
                    h0=use_h0,four_char=four_char,cosmological=cosmological,skip_bh=skip_bh, \
                    GAS=SHOW_GAS,BOX_MAX=xmax_tmp,CENTER_BOX=center_i, \
                    min_stellar_age=min_stellar_age);
            if (SHOW_STARS==1):
                id_i_g,m_i_g,x_i_g,y_i_g,z_i_g,vx_i_g,vy_i_g,vz_i_g,h_i_g,c_i_g,zm_i_g = \
                  load_pos_vel_etc(snapdir+'/'+subfolder_ext, snapshot_list[i_snap], filename_prefix = snapdir_specific, \
                    h0=use_h0,four_char=four_char,cosmological=cosmological,skip_bh=skip_bh, \
                    GAS=1,SWAP_TEMP_RHO=1,BOX_MAX=xmax_tmp,CENTER_BOX=center_i, \
                    min_stellar_age=min_stellar_age);
            ## correct to comoving coordinates for interpolation (so correctly capture hubble flow)
            ## (put this here, splitting initial and final rescalings so dont do it twice on 'recycling' step)
            if (cosmological==1):
                for vec in [x_i,y_i,z_i,h_i,vx_i,vy_i,vz_i]: vec /= a_scale_i;
                if (SHOW_STARS==1):
                    for vec in [x_i_g,y_i_g,z_i_g,h_i_g,vx_i_g,vy_i_g,vz_i_g]: vec /= a_scale_i;

            ## now load snapshot (f) [should be new, have to load fresh] 
            id_f,m_f,x_f,y_f,z_f,vx_f,vy_f,vz_f,h_f,c_f,zm_f = \
                load_pos_vel_etc(snapdir+'/'+subfolder_ext, snapshot_list[i_snap+1], filename_prefix = snapdir_specific, \
                    h0=use_h0,four_char=four_char,cosmological=cosmological,skip_bh=skip_bh, \
                    GAS=SHOW_GAS,BOX_MAX=xmax_tmp,CENTER_BOX=center_f);
            if (SHOW_STARS==1):
                id_f_g,m_f_g,x_f_g,y_f_g,z_f_g,vx_f_g,vy_f_g,vz_f_g,h_f_g,c_f_g,zm_f_g = \
                  load_pos_vel_etc(snapdir+'/'+subfolder_ext, snapshot_list[i_snap+1], filename_prefix = snapdir_specific, \
                    h0=use_h0,four_char=four_char,cosmological=cosmological,skip_bh=skip_bh, \
                    GAS=1,SWAP_TEMP_RHO=1,BOX_MAX=xmax_tmp,CENTER_BOX=center_f);
            ## correct to comoving coordinates for interpolation (so correctly capture hubble flow)
            if (cosmological==1):
                for vec in [x_f,y_f,z_f,h_f,vx_f,vy_f,vz_f]: vec /= a_scale_f;
                if (SHOW_STARS==1):
                    for vec in [x_f_g,y_f_g,z_f_g,h_f_g,vx_f_g,vy_f_g,vz_f_g]: vec /= a_scale_f;
            snapshot_f_number_that_i_just_loaded = snapshot_list[i_snap+1];

            ## before going further, check that we actually have particles to process, 
            ##   and if not skip this frame
            if (m_i.size+m_f.size < 5):
                continue;
            
            
            print('... ... matching ids in the snapshots ... ...')
            ## predefine matches to save time in the loop below
            sub_id_i,sub_id_f,nomatch_from_i_in_f,nomatch_from_f_in_i = compile_matched_ids( id_i, id_f );
            if (SHOW_STARS==1):
                sub_id_i_g,sub_id_f_g,nomatch_from_i_in_f_g,nomatch_from_f_in_i_g = compile_matched_ids( id_i_g, id_f_g );
    
            print('... ... entering frame loop ... ...')
            ## loop over frames in between the two snapshots just pulled up
            for i_frame in range(n_frames_to_do):
                j_of_frame = frames_to_do[i_frame];
                dt = time_frame_grid[j_of_frame] - t_i
                time_of_frame = time_frame_grid[j_of_frame]

                print('... ... ... interpolating for frame ',i_frame+1,'/',n_frames_to_do,'... ... ...')
                ## this is the interpolation step between the two snapshots for each frame
                x_all,y_all,z_all,m_all,h_all,c_all,zm_all = interpolation_for_movie( dt, delta_t, \
                    m_i,h_i,c_i,zm_i, m_f,h_f,c_f,zm_f, \
                    x_i,y_i,z_i,vx_i,vy_i,vz_i, x_f,y_f,z_f,vx_f,vy_f,vz_f, \
                    sub_id_i,sub_id_f, nomatch_from_i_in_f,nomatch_from_f_in_i, \
                    center_i,center_f,use_polar_interpolation=use_polar_interpolation, set_dx_frame=2.*x_plot_scale_grid[j_of_frame], stars=SHOW_STARS )
                if (SHOW_STARS==1):
                    x_all_g,y_all_g,z_all_g,m_all_g,h_all_g,c_all_g,zm_all_g = \
                        interpolation_for_movie( dt, delta_t, \
                        m_i_g,h_i_g,c_i_g,zm_i_g,m_f_g,h_f_g,c_f_g,zm_f_g, \
                        x_i_g,y_i_g,z_i_g,vx_i_g,vy_i_g,vz_i_g, x_f_g,y_f_g,z_f_g,vx_f_g,vy_f_g,vz_f_g, \
                        sub_id_i_g,sub_id_f_g, nomatch_from_i_in_f_g,nomatch_from_f_in_i_g, \
                        center_i,center_f,use_polar_interpolation=use_polar_interpolation, set_dx_frame=2.*x_plot_scale_grid[j_of_frame], stars=0 )
                ## ok, now have the interpolated gas(+stellar) quantities needed for each image

                ## correct back from comoving coordinates for plotting (since units are physical)
                if (cosmological==1):
                    a = a_scale_grid[j_of_frame]
                    time_of_frame = a
                    for vec in [x_all,y_all,z_all,h_all]: vec *= a;
                    if (SHOW_STARS==1):
                        for vec in [x_all_g,y_all_g,z_all_g,h_all_g]: vec *= a;

                print('... ... ... centering and setting passed variables ... ... ...')
                ## re-center on the current frame center, interpolated between the snapshots
                cen = get_precalc_zoom_center(snapdir+'/'+subfolder_ext, snapdir_specific, outputdir_master,a_scale_grid[j_of_frame],cosmological=cosmological)
                print('Center for frame j=',j_of_frame,' at ',cen)
                x_all-=cen[0]; y_all-=cen[1]; z_all-=cen[2];
                if (SHOW_STARS==1):
                    x_all_g-=cen[0]; y_all_g-=cen[1]; z_all_g-=cen[2];
        
    
                ## set dynamic ranges of image
                maxden=maxden_grid[j_of_frame]
                dynrange=dynrange_grid[j_of_frame]
                ## set image spatial scale 
                xr_0=x_plot_scale_grid[j_of_frame]
                yr_0=y_plot_scale_grid[j_of_frame]
                xr=[-xr_0,xr_0]; yr=[-yr_0,yr_0]; zr_0=max([xr_0,yr_0]); zr=[-zr_0,zr_0]
                ## set viewing angle for image
                theta = theta_frame_grid[j_of_frame]
                phi   = phi_frame_grid[j_of_frame]

                ## some final variable cleaning/prep before imaging routine:
                if (SHOW_STARS==0):
                    gdum=np.array([0]);
                    x_all_g=y_all_g=z_all_g=m_all_g=h_all_g=c_all_g=zm_all_g=gdum;
                else:
                    gdum=np.zeros(m_all_g.size)+1.;
                        
                print('... ... ... sending to main imaging routine ... ... ...')
                if(2==0):
                    print(snapdir_specific)
                    print(j_of_frame)
                    print(snapdir_master)
                    print(outputdir_master)
                    print(theta)
                    print(phi)
                    print(dynrange)
                    print(maxden)
                    print(show_gasstarxray)
                    print(add_extension)
                    print(show_time_label)
                    print(show_scale_label)
                    print(use_h0)
                    print(cosmological)
                    print(pixels)
                    print(xr)
                    print(yr)
                    print(zr)
                    print(nasa_colors)
                    print(sdss_colors)
                    print(use_old_extinction_routine)
                    print(do_with_colors)
                    print(threecolor)
                    print(temp_min)
                    print(temp_max)
                    print(np.array([temp_min, temp_max]))
                    print(SHOW_GAS)
                    print(m_all)
                    print(x_all)
                    print(y_all)
                    print(z_all)
                    print(c_all)
                    print(h_all)
                    print(zm_all)
                    print(gdum)
                    print(c_all_g)
                    print(h_all_g)
                    print(zm_all_g)
                    print(m_all_g)
                    print(x_all_g)
                    print(y_all_g)
                    print(z_all_g)
                    print(time_of_frame)

                ## alright, now we can actually send this to the imaging routine 
                ##   to make the frame image! 
                image24, massmap = \
                  image_maker( snapdir_specific, j_of_frame, \
                    snapdir_master=snapdir_master, outdir_master=outputdir_master, 
                    theta=theta, phi=phi, dynrange=dynrange, maxden=maxden, 
                    show_gasstarxray=show_gasstarxray, add_extension=add_extension, \
                    show_time_label=show_time_label, show_scale_label=show_scale_label, \
                    use_h0=use_h0, cosmo=cosmological, coordinates_cylindrical=0, \
                    set_percent_maxden=0, set_percent_minden=0, pixels=pixels, \
                    center_on_com=0, center_on_bh=0, center=[0., 0., 0.], \
                    xrange=xr, yrange=yr, zrange=zr, \
                    project_to_camera=1, camera_opening_angle=camera_opening_angle, \
                    center_is_camera_position=0, camera_direction=[0.,0.,-1.], \
                    nasa_colors=nasa_colors, sdss_colors=sdss_colors, \
                    use_old_extinction_routine=use_old_extinction_routine, \
                    do_with_colors=do_with_colors, threecolor=threecolor, \
                    log_temp_wt=1, set_temp_min=temp_min, set_temp_max=temp_max, \
                    gas_map_temperature_cuts=np.array([temp_min, temp_max]), \
                    input_data_is_sent_directly=1, include_lighting=include_lighting, \
                    m_all=m_all,x_all=x_all,y_all=y_all,z_all=z_all,c_all=c_all,h_all=h_all,zm_all=zm_all+z_to_add,\
                    gas_u=gdum,gas_rho=c_all_g,gas_numh=gdum,gas_nume=gdum,gas_hsml=h_all_g,gas_metallicity=zm_all_g+z_to_add,\
                    gas_mass=m_all_g,gas_x=x_all_g,gas_y=y_all_g,gas_z=z_all_g,time=time_of_frame, \
                    scattered_fraction=scattered_fraction,min_stellar_age=min_stellar_age, \
                    snapshot_subdirectory=subfolder_ext,set_added_layer_ctable=set_added_layer_ctable, \
                    add_gas_layer_to_image=add_gas_layer_to_image, \
                    dust_to_gas_ratio_rescale=dust_to_gas_ratio_rescale,set_added_layer_alpha=set_added_layer_alpha, \
                    invert_colors=invert_colors,spin_movie=spin_movie, \
                    h_rescale_factor=h_rescale_factor,h_max=h_max, \
                    projection_vector_x=projection_vector_x,projection_vector_y=projection_vector_y,projection_vector_z=projection_vector_z);        


    

                
                ## don't need to do anything else here, the image & datafile are dumped!
                print('... ... ... frame ',j_of_frame,'complete! ... ... ...')
              
    gc.collect()
    return 1; ## success!


def image_maker( sdir, snapnum, \
    snapdir_master='/n/scratch2/hernquist_lab/phopkins/sbw_tests/', 
    outdir_master='/n/scratch2/hernquist_lab/phopkins/images/', 
    snapshot_subdirectory='output',
    theta=0., phi=0., dynrange=1.0e5, maxden_rescale=1., maxden=0.,
    show_gasstarxray = 'gas', #	or 'star' or 'xray'
    add_gas_layer_to_image='', set_added_layer_alpha=0.3, set_added_layer_ctable='heat_purple', 
    add_extension='', show_time_label=1, show_scale_label=1, 
    filename_set_manually='',
    do_with_colors=1, log_temp_wt=1, include_lighting=1, 
    set_percent_maxden=0, set_percent_minden=0, 
    center_on_com=0, center_on_bh=0, use_h0=1, cosmo=0, 
    center=[0., 0., 0.], coordinates_cylindrical=0, 
    offset=[0.,0.,0.],
    pixels=720,xrange=[-1.,1.],yrange=0,zrange=0,set_temp_max=0,set_temp_min=0,
    set_size_inches=1.,axis_enabled=False,
    threecolor=1, nasa_colors=1, sdss_colors=0, use_old_extinction_routine=0, 
    project_to_camera=1, camera_opening_angle=45.0,
    center_is_camera_position=0, camera_direction=[0.,0.,-1.], 
    gas_map_temperature_cuts=[1.0e4, 1.0e6], 
    input_data_is_sent_directly=0,
    m_all=0,x_all=0,y_all=0,z_all=0,c_all=0,h_all=0,zm_all=0,
    gas_u=0,gas_rho=0,gas_hsml=0,gas_numh=0,gas_nume=0,gas_metallicity=0,
    gas_mass=0,gas_x=0,gas_y=0,gas_z=0,time=0,
    invert_colors=0,spin_movie=0,scattered_fraction=0.01,
    min_stellar_age=0,h_rescale_factor=1.,h_max=0,
    angle_to_rotate_image=0.,
    projection_vector_x=[1,0,0],projection_vector_y=[0,1,0],projection_vector_z=[0,0,1],
    BAND_IDS=[9,10,11], ## ugr composite
    set_f_saturated=0.0001,
    dust_to_gas_ratio_rescale=1.0):
    

	## define some of the variables to be used
    ss=snap_ext(snapnum,four_char=1);
    tt=snap_ext(np.around(theta).astype(int));
    theta *= np.pi/180.; phi *= np.pi/180.; # to radians
    nameroot = sdir+'_s'+ss+'_t'+tt;
    if(spin_movie==1):
        nameroot = sdir+'_s'+ss+'_spin';    
    outputdir = outdir_master+sdir;
    call(["mkdir",outputdir]);
    outputdir+='/'; 
    snapdir=snapdir_master+sdir+'/'+snapshot_subdirectory+'/';
    suff='_'+show_gasstarxray;
    if (threecolor==1):
        if (sdss_colors==1): suff+='_S3c'
        if (nasa_colors==1): suff+='_N3c'
    suff+=add_extension;
    fname_base=outputdir+nameroot+suff
    do_xray=0; do_stars=0; 
    if((show_gasstarxray=='xr') or (show_gasstarxray=='xray')): do_xray=1; do_with_colors=0;
    if((show_gasstarxray=='star') or (show_gasstarxray=='stars') or (show_gasstarxray=='st')): do_stars=1;
    
    
    ## check whether the data needs to be pulled up, or if its being passed by the calling routine
    if(input_data_is_sent_directly==0):
        ## read in snapshot data and do centering 
        ptypes=[2,3,4]; ## stars in non-cosmological snapshot
        if (cosmo==1): ptypes=[4];
        if (do_stars==0): ptypes=[0];
        m_all, x_all, y_all, z_all, c_all, h_all, zm_all, time = load_snapshot_brick(ptypes, \
            snapdir, snapnum, h0=use_h0, cosmological=cosmo, filename_prefix=sdir, \
            skip_bh=cosmo, do_xray=do_xray, four_char=0, use_rundir=0, full_gas_set=0);
        h_all *= 1.25;

        if ((do_stars==1) & (threecolor==1)): ## will need gas info to process attenuation
            gas_u, gas_rho, gas_hsml, gas_numh, gas_nume, gas_sfr, gas_lxray, gas_metallicity, gas_mass, gas_x, gas_y, gas_z, time = \
             load_snapshot_brick([0], snapdir, snapnum, h0=use_h0, cosmological=cosmo, full_gas_set=1);
            ## don't allow hot gas to have dust
            gas_temp = gadget.gas_temperature(gas_u,gas_nume); gas_metallicity[gas_temp > 1.0e6] = 0.0; 
    
        if (center[0] == 0.):
            if (cosmo==1): 
                center,n00 = fast_zoom_center(snapdir,snapnum,cosmological=cosmo);
            else:
                if (center_on_bh==1):
                    pbh=gadget.readsnap(snapdir,snapnum,5,h0=use_h0,cosmological=0,skip_bh=0);
                    pos=pbh['p']; center=[pos[0,0],pos[0,1],pos[0,2]];
                if (center_on_com==1):
                    center=[np.median(x_all),np.median(y_all),np.median(z_all)];
        center=np.array(center)+np.array(offset)
        x_all-=center[0]; y_all-=center[1]; z_all-=center[2];
        print('center at ',center)


    ## rotate and re-project coordinate frame 
    if (center_is_camera_position==0):
        x=1.*x_all; y=1.*y_all; z=1.*z_all; 
        #x,y,z=coordinates_rotate(x_all,y_all,z_all,theta,phi,coordinates_cylindrical=coordinates_cylindrical);
        if ((do_stars==1) & (threecolor==1)):
            gas_x-=center[0]; gas_y-=center[1]; gas_z-=center[2];
        gx=gas_x; gy=gas_y; gz=gas_z;
            #gx,gy,gz=coordinates_rotate(gas_x,gas_y,gas_z,theta,phi,coordinates_cylindrical=coordinates_cylindrical);
    else:
        x=x_all; y=y_all; z=z_all; 
        if ((do_stars==1) & (threecolor==1)): gx=gas_x; gy=gas_y; gz=gas_z;
    
    x_tmp=1.*x; y_tmp=1.*y; z_tmp=1.*z;
    x=projection_vector_x[0]*x_tmp + projection_vector_x[1]*y_tmp + projection_vector_x[2]*z_tmp
    y=projection_vector_y[0]*x_tmp + projection_vector_y[1]*y_tmp + projection_vector_y[2]*z_tmp
    z=projection_vector_z[0]*x_tmp + projection_vector_z[1]*y_tmp + projection_vector_z[2]*z_tmp
    if ((do_stars==1) & (threecolor==1)): 
        x_tmp=1.*gx; y_tmp=1.*gy; z_tmp=1.*gz;
        gx=projection_vector_x[0]*x_tmp + projection_vector_x[1]*y_tmp + projection_vector_x[2]*z_tmp
        gy=projection_vector_y[0]*x_tmp + projection_vector_y[1]*y_tmp + projection_vector_y[2]*z_tmp
        gz=projection_vector_z[0]*x_tmp + projection_vector_z[1]*y_tmp + projection_vector_z[2]*z_tmp

    if(center_is_camera_position==0):
        x,y,z=coordinates_rotate(x,y,z,theta,phi,coordinates_cylindrical=coordinates_cylindrical);
        if ((do_stars==1) & (threecolor==1)):
            gx,gy,gz=coordinates_rotate(gx,gy,gz,theta,phi,coordinates_cylindrical=coordinates_cylindrical);



    ## set dynamic ranges of image
    temp_max=1.0e6; temp_min=1.0e3; ## max/min gas temperature
    if (do_stars==1): temp_max=50.; temp_min=0.01; ## max/min stellar age
    if (set_temp_max != 0): temp_max=set_temp_max
    if (set_temp_min != 0): temp_min=set_temp_min	
    xr=xrange; yr=xr; zr=0;
    if (checklen(yrange)>1): yr=yrange;
    if (checklen(zrange)>1): zr=zrange;
    scale=0.5*np.max([xr[1]-xr[0],yr[1]-yr[0]]); zr=np.array([-1.,1.])*scale;
    if(maxden==0.):
        maxden = 0.1 * 1.e11/1.0e10 * (2.0/scale)**(0.3) * maxden_rescale;
    if ((threecolor==1) & (do_stars)):
        if(dynrange==1.0e5): # reset if it's the default value 
            dynrange=1.0e2;
            if (nasa_colors==1): dynrange=1.0e4;
            if (nasa_colors==1): maxden *= 30.;


    ## now check if we're projecting to a camera (instead of a fixed-plane)
    xlen=0.5*(xr[1]-xr[0]); ylen=0.5*(yr[1]-yr[0]); 
    xr_0=xr; yr_0=yr;
    if (project_to_camera==1):  
        ## determine the angular opening and camera distance given physical x/y range:
        ## note: camera_opening_angle is the half-angle
        xr=np.array([-1.,1.])*np.tan(camera_opening_angle*np.pi/180.); yr=xr*ylen/xlen; ## n- degree opening angle
        ## use this to position camera
        c_dist = xlen/xr[1]; 
        camera_direction /= np.sqrt(camera_direction[0]**2.+camera_direction[1]**2.+camera_direction[2]**2.);
        camera_pos = -c_dist * camera_direction;
        ## and determine z-depth of image
        cpmax=np.sqrt(camera_pos[0]**2.+camera_pos[1]**2.+camera_pos[2]**2.); zrm=np.sqrt(zr[0]**2.+zr[1]**2.);
        ## clip z_range for speed : 
        zrm=5.0*zrm; 
        #zrm=2.5*zrm
        zrm=np.sqrt(zrm*zrm+cpmax*cpmax); zr=[-zrm,0.]
        ## correct if the center given is the camera position:
        if (center_is_camera_position==1): camera_pos=[0.,0.,0.];

        ## now do the actual projection into camera coordinates:        
        #r_pre_camera = np.sqrt(x*x+y*y+z*z)
        x,y,z,rc=coordinates_project_to_camera(x,y,z,camera_pos=camera_pos,camera_dir=camera_direction);
        h_all*=rc; ## correct for size re-scaling

        #TEMPORARY_FOR_SPIN_MOVIE_ONLY = spin_movie
        #if(TEMPORARY_FOR_SPIN_MOVIE_ONLY==1):
        #    h_max_camera = 0.1 * 0.01 * (1. + (np.abs(z/(0.5*c_dist)))**5.) # too 'grey'? (copy 3, nothing else)
        #    h_max_camera = np.maximum(0.005+0.*z , 0.1 * 0.01 * (1. + (np.abs(z/(0.5*c_dist)))**5.)) # copy 2 
        #    h_max_camera = np.maximum(0.00707+0.*z , 0.1 * 0.01 * (1. + (np.abs(z/(0.5*c_dist)))**5.)) # copy 6
        #    h_max_camera = 0.01 * (1. + (np.abs(z/(0.5*c_dist)))**5.) # good foreground; copy 7
            #h_max_camera = np.minimum(0.01 * (1. + (np.abs(z/(0.5*c_dist)))**5.), 0.01 * (1. + (np.abs(z/(0.5*c_dist)))**0.5)) # copy 4
            #h_max_camera *= 0.7 * np.minimum((1.+(c_all/2.)**2.),(1.+(c_all/5.0)**0.75)) # copy 5 (with copy 4 line also)
            #h_max_camera *= 1. + (r_pre_camera/10.)**0.5
            #h_max_camera = np.maximum(0.005*(1.+3.0/np.abs(z/(0.1*c_dist))) , 0.1 * 0.01 * (1. + (np.abs(z/(0.5*c_dist)))**5.)) # copy 8 (1of2)
            #h_max_camera = np.maximum(0.0071*(1.+1.5/np.abs(z/(0.1*c_dist))**2.) , 0.1 * 0.01 * (1. + (np.abs(z/(0.5*c_dist)))**5.)) # copy 9 (1of2)
            #h_max_camera = np.maximum(0.0071*(1.+1.5/np.abs(z/(0.1*c_dist))**2.) , 0.1 * 0.01 * (1. + (np.abs(z/(0.5*c_dist)))**2.)) # copy 9 good except too 'sharp' a transition to smooth for later snaps (b/c of ~5 power); lower here (copy 10)
            #h_max_camera = np.maximum(0.0071/((1.+scale/2.)**0.125)*(1.+1.5/np.abs(z/(0.1*c_dist))**2.) , 0.1 * 0.01 * (1. + (np.abs(z/(0.5*c_dist)))**2.)) # copy 9 good except too 'sharp' a transition to smooth for later snaps (b/c of ~5 power); lower here (copy 11)
            #h_max_camera = (0.0071/((1.+scale/2.)**0.125)*(0.8+0.7/np.abs(z/(0.1*c_dist))**2.)) * (1. + (np.abs(z/(0.5*c_dist)))*2.5) # good foreground; copy 7
        #    h_max_camera = 0.01 * (1. + (np.abs(z/(0.5*c_dist)))**2.) # good foreground; copy 7_pt5 set
        #    h_all = np.minimum(h_all,h_max_camera)
        #    m_all*=1./(z*z + h_all*h_all + (0.*0.25*c_dist)**2.); ## makes into 'fluxes' (dimmer when further away)
            #m_all *= (1. + ((0.02/(0.0001+h_all))**2.) / (1. + np.abs(z/(0.1*c_dist))**2. ));  # copy 8 (2of2)
            #m_all *= (1. + ((0.02/(0.0001+h_all))**2.) / (1. + np.abs(z/(0.1*c_dist))**2. ));  # copy 9/10/11 (2of2)
            #m_all *= (1. + ((0.02/(0.0001+h_all))**2.) / (1. + np.abs(z/(0.1*c_dist))**2. ));  # copy 9/10/11 (2of2)
        #else:

        #h_max_camera = 0.02 * (h_all/0.02)**(0.5)
        #h_max_camera = 0.01 * (1. + (np.abs(z/(0.5*c_dist)))**(-2.)) # good foreground; copy 7_pt5 set
        #h_all = h_all * (h_all / 0.01)**(0.5)
        h_all_lim = 0.01
        biggie=np.where(h_all > h_all_lim); h_all[biggie] *= (h_all[biggie]/h_all_lim)**(0.0)
        #h_all = h_all * (h_all / 0.01)**(-0.5)
        h_max_camera = 1.*h_all
        #biggie = np.where(h_all > h_max_camera)
        #m_all[biggie] *= (h_max_camera[biggie]/h_all[biggie])**(0.5)
        h_all = np.minimum(h_all,h_max_camera)
        

        m_all*=1./(z*z + h_all*h_all + (0.25*c_dist)**2.); ## makes into 'fluxes' (dimmer when further away)
        ## clip particles too close to camera to prevent artifacts
        z_clip = -c_dist/10.0
        m_all[z >= z_clip] = 0.0
        ## also need to correct gas properties if using it for attenuation
        if ((threecolor==1) & (do_stars==1)):
            gx,gy,gz,grc=coordinates_project_to_camera(gx,gy,gz,camera_pos=camera_pos,camera_dir=camera_direction);
            gas_hsml*=grc; gas_mass*=1./(gz*gz + gas_hsml*gas_hsml + (0.25*c_dist)**2.); ## this conserves surface density, so attenuation units are fine
            gas_mass[gz >= z_clip] = 0.0 ## clipping to prevent artifacts
            biggie=np.where(gas_hsml > h_all_lim); gas_hsml[biggie] *= (gas_hsml[biggie]/h_all_lim)**(0.0)
        
        #rotate the image; assumes angle is specified in degrees
        angle_to_rotate_image*=np.pi/180.
        xnew = x*np.cos(angle_to_rotate_image) - y*np.sin(angle_to_rotate_image)
        ynew = x*np.sin(angle_to_rotate_image) + y*np.cos(angle_to_rotate_image)
        x = xnew
        y = ynew
        if ((threecolor==1) & (do_stars==1)):
            gxnew = gx*np.cos(angle_to_rotate_image) - gy*np.sin(angle_to_rotate_image)
            gynew = gx*np.sin(angle_to_rotate_image) + gy*np.cos(angle_to_rotate_image)
            gx = gxnew
            gy = gynew

    #plotting_setup_junk
    plt.close('all')
    format = '.pdf' #format = '.png' ## '.ps','.pdf','.png' also work well
    axis_ratio = ylen/xlen
    #fig=plt.figure(frameon=False,figsize=(10.,10.*axis_ratio),dpi=100)
    ## for Chris's high-res images:
    fig=plt.figure(frameon=False,dpi=pixels,figsize=(1.,1.*axis_ratio))
    fig.set_size_inches(set_size_inches,set_size_inches*axis_ratio)
    ax_fig=plt.Axes(fig,[0.,0.,1.,1.*axis_ratio],clip_on=True)
    if(axis_enabled==False): ax_fig.set_axis_off()
    frame1=plt.gca()
    if(axis_enabled==False): frame1.set_axis_off()
    plt.xlabel(''); plt.ylabel('')
    frame1.axes.xaxis.set_ticklabels([]) #no tick names
    frame1.axes.yaxis.set_ticklabels([]) #no tick names
    frame1.axes.get_xaxis().set_ticks([]) #no ticks
    frame1.axes.get_yaxis().set_ticks([]) #no ticks
    ax_fig.axes.xaxis.set_ticklabels([]) #no tick names
    ax_fig.axes.yaxis.set_ticklabels([]) #no tick names
    ax_fig.axes.get_xaxis().set_ticks([]) #no ticks
    ax_fig.axes.get_yaxis().set_ticks([]) #no ticks
    if(axis_enabled==True): fig.add_axes(ax_fig)
    if(axis_enabled==False): ax_fig=frame1

	  
	## trim particles to those inside/near the plotting region
    sidebuf=10.; 
    if (cosmo==1): sidebuf=1.15
    x=np.array(x,dtype='f'); y=np.array(y,dtype='f'); z=np.array(z,dtype='f'); 
    m_all=np.array(m_all,dtype='f'); h_all=np.array(h_all,dtype='f'); 
    c_all=np.array(c_all,dtype='f'); 
    ok =    ok_scan(x,xmax=np.max(np.fabs(xr))*sidebuf) & ok_scan(y,xmax=np.max(np.fabs(yr))*sidebuf) & \
            ok_scan(z,xmax=np.max(np.fabs(zr))*sidebuf) & ok_scan(m_all,pos=1,xmax=1.0e40) & \
            ok_scan(h_all,pos=1) & ok_scan(c_all,pos=1,xmax=1.0e40)
    weights = m_all; color_weights = c_all; ## recall, gadget masses in units of 1.0d10
    print("Mean h_all = ",h_all.mean())
    print("Median h_all = ",np.median(h_all))


    TEMPORARY_FOR_THORSTEN_REVIEW_ONLY = 0
    if(TEMPORARY_FOR_THORSTEN_REVIEW_ONLY==1):
        #xxx=1.*x; x=1.*y; y=1.*xxx;
        if ((threecolor==1) & (do_stars==1)):
            #gxxx=1.*gx; gx=1.*gy; gy=1.*gxxx;
            h_all *= 0.8
            #c_all *= (c_all / 1.0);
            gas_metallicity *= 1.e-5;
            c_all *= (c_all / 1.0)**0.5; gas_metallicity *= 0.67;
            #c_all *= (c_all / 1.0)**0.25; gas_metallicity *= 0.67;
            #c_all *= (c_all / 1.0)**0.5; gas_metallicity *= 0.67;
            #c_all *= (c_all / 1.0)**0.25; 
            use_old_extinction_routine = 0;

    TEMPORARY_FOR_SPIN_MOVIE_ONLY = spin_movie
    if(TEMPORARY_FOR_SPIN_MOVIE_ONLY==1):
        if ((threecolor==1) & (do_stars==1)):
            h_all = np.maximum(h_all, 1.e-40)
            h_all *= ((scale/30.)**0.25) * (1.0 + (h_all / 0.1)**(0.5)) / (1. + (h_all / 1.0)**(0.5))
            #c_all *= (c_all / 1.0)**0.25; 
            c_all *= (c_all / 1.0)**0.125; 
            #h_all *= h_rescale_factor
            ## limit softenings to a maximum fraction of FOV to avoid foreground fuzziness/blobs
            #h_all = np.minimum(h_all,h_max) # CCH
        else:
            h_all *= (((30.*0.+1.*scale)/30.)**0.25) * (0.8 + (h_all / 0.1)**(0.5)) / (1. + (h_all / 1.0)**(0.5))

    ## adjust softenings if appropriate keywords are set
    h_all *= h_rescale_factor
    if (h_max > 0):
        h_all = np.minimum(h_all,h_max)

    ## alright, now ready for the main plot construction:
    if (threecolor==1):
        if (do_stars==1):
            ## making a mock 3-color composite image here ::
            
            ## first grab the appropriate, attenuated stellar luminosities
            ## for now not including BHs, so initialize some appropriate dummies:
            bh_pos=[0,0,0]; bh_luminosity=0.; include_bh=0;
            ## and some limits for the sake of integration:
            SizeMax=20.*np.max([np.fabs(xr[1]-xr[0]),np.fabs(yr[1]-yr[0])]);
            SizeMin=np.min([np.fabs(xr[1]-xr[0]),np.fabs(yr[1]-yr[0])])/250.;
 
            ## now set up the main variables:
            star_pos=np.zeros((3,checklen(x[ok]))); 
            star_pos[0,:]=x[ok]; star_pos[1,:]=y[ok]; star_pos[2,:]=z[ok];
            stellar_age=np.maximum(c_all,min_stellar_age);
            stellar_metallicity=zm_all; stellar_mass=m_all;
            print('Min stellar age = ',stellar_age.min())
            print('Mean stellar age = ',stellar_age.mean())
            print('Mean stellar metallicity = ',stellar_metallicity.mean())
            ## gas will get clipped in attenuation routine, don't need to worry about it here
            gas_pos=np.zeros((3,checklen(gx))); 
            gas_pos[0,:]=gx; gas_pos[1,:]=gy; gas_pos[2,:]=gz;

            if ((add_gas_layer_to_image=='')==False):
                gas_wt = 0.*gx
                maxden_for_layer=1.*maxden; dynrange_for_layer=1.*dynrange
                set_percent_maxden_layer=0.; set_percent_minden_layer=0.;
                gas_hsml_for_extra_layer=gas_hsml
                if (add_gas_layer_to_image=='Halpha'): 
                    gas_wt = gas_mass * gas_rho*gas_nume*(1.-gas_numh)*(gadget.gas_temperature(gas_u,gas_nume)**(-0.75))/(gadget.gas_mu(gas_nume)**2.)
                    maxden_for_layer *= 3.e4;
                    #dynrange_for_layer *= 1.5;
                    dynrange_for_layer *= 0.7;
                    dynrange_for_layer *= 2.0;
                    #gas_hsml_for_extra_layer *= 1.1;
                if (add_gas_layer_to_image=='CO'): 
                    gas_tmp = gadget.gas_temperature(gas_u,gas_nume)
                    n_tmp = gas_rho * 176.2; # assuming mean molec weight of 2.3 for dense gas
                    gas_wt = gas_mass * np.exp(-(gas_tmp/8000. + 10./n_tmp));
                    maxden_for_layer *= 3.e4;
                    dynrange_for_layer *= 0.7e4;
                if (add_gas_layer_to_image=='SFR'): 
                    gas_wt = gas_sfr
                    set_percent_maxden_layer=0.9999;
                    set_percent_minden_layer=0.01;
                if (add_gas_layer_to_image=='Xray'): 
                    gas_wt = gas_lxray
                    maxden_for_layer *= 0.6e3;
                    dynrange_for_layer *= 0.1;
                    gas_hsml_for_extra_layer *= 1.2;
                if (add_gas_layer_to_image=='Zmetal'): 
                    gas_wt = gas_mass*gas_metallicity
                    set_percent_maxden_layer=0.9999;
                    set_percent_minden_layer=0.01;
                gas_wt /= np.sum(gas_wt)
                                    
                massmap_gas_extra_layer,image_singledepth_extra_layer = \
                cmakepic.simple_makepic(gx,gy,weights=gas_wt,hsml=gas_hsml_for_extra_layer,\
                    xrange=xr,yrange=yr,
                    set_dynrng=dynrange_for_layer,set_maxden=maxden_for_layer,
                    set_percent_maxden=set_percent_maxden_layer,set_percent_minden=set_percent_minden_layer, 
                    color_temperature=0,pixels=pixels,invert_colorscale=1-invert_colors);


            if (use_old_extinction_routine == 0):
                ##
                ## this is the newer 'single call' attenuation and ray-tracing package: 
                ##  slightly more accurate in diffuse regions, more stable behavior
                ##

                if (x[ok].size <= 3):
                    ## we got a bad return from the previous routine, initialize blank arrays
                    out_gas=out_u=out_g=out_r=np.zeros((pixels,pixels))
                    image24=massmap=np.zeros((pixels,pixels,3))
                else:
                    ## actually call the ray-trace:
                    out_gas,out_u,out_g,out_r = rayproj.stellar_raytrace( BAND_IDS, \
                        x[ok], y[ok], z[ok], \
                        stellar_mass[ok], stellar_age[ok], stellar_metallicity[ok], h_all[ok], \
                        gx, gy, gz, gas_mass, \
                        gas_metallicity*dust_to_gas_ratio_rescale, gas_hsml, \
                        xrange=xr, yrange=yr, zrange=zr, pixels=pixels, \
                        #ADD_BASE_METALLICITY=1.0*0.02, ADD_BASE_AGE=0.003, 
                        ADD_BASE_METALLICITY=1.0e-6*0.02, ADD_BASE_AGE=0.0003, 
                        #ADD_BASE_METALLICITY=0.01*0.02, ADD_BASE_AGE=0.0003, 
                        IMF_SALPETER=0, IMF_CHABRIER=1 );
                
                    if(np.array(out_gas).size<=1):
                        ## we got a bad return from the previous routine, initialize blank arrays
                        out_gas=out_u=out_g=out_r=np.zeros((pixels,pixels))
                        image24=massmap=np.zeros((pixels,pixels,3))
                    else:
                        ## make the resulting maps into an image
                        image24, massmap, image_max, image_min = \
                            makethreepic.make_threeband_image_process_bandmaps( out_r,out_g,out_u, \
                            maxden=maxden,dynrange=dynrange,pixels=pixels, \
                            color_scheme_nasa=nasa_colors,color_scheme_sdss=sdss_colors, \
                            set_f_saturated=set_f_saturated);  
                
            else: ## use_old_extinction_routine==1; can be useful for certain types of images
                ## (for example, this can allow for a scattered fraction, the other does not currently)

                ## call routine to get the post-attenuation band-specific luminosities
                lum_noatten, lum_losNH = \
                getlum.get_attenuated_stellar_luminosities( BAND_IDS, star_pos, gas_pos, bh_pos, \
                    stellar_age[ok], stellar_metallicity[ok], stellar_mass[ok], \
                    gas_u, gas_rho, gas_hsml, gas_numh, gas_nume, \
                    gas_metallicity*dust_to_gas_ratio_rescale, gas_mass, \
                    bh_luminosity, \
                    xrange=xr, yrange=yr, zrange=zr, \
                    INCLUDE_BH=include_bh, SKIP_ATTENUATION=0, 
                    #ADD_BASE_METALLICITY=1.0, ADD_BASE_AGE=0., 
                    ADD_BASE_METALLICITY=1.0e-2, ADD_BASE_AGE=0., 
                    IMF_SALPETER=0, IMF_CHABRIER=1, \
                    MIN_CELL_SIZE=SizeMin, OUTER_RANGE_OF_INT=SizeMax, \
                    SCATTERED_FRACTION=scattered_fraction, \
                    REDDENING_SMC=1, REDDENING_LMC=0, REDDENING_MW=0, \
                    AGN_MARCONI=0, AGN_HRH=1, AGN_RICHARDS=0, AGN_SDSS=0 );
                
                ## now pass these into the 3-band image routine
                image24, massmap = \
                makethreepic.make_threeband_image(x[ok],y[ok],lum_losNH,hsml=h_all[ok],\
                    xrange=xr,yrange=yr,maxden=maxden,dynrange=dynrange,pixels=pixels, \
                    color_scheme_nasa=nasa_colors,color_scheme_sdss=sdss_colors );           

        else: 
            ## 
            ## threecolor==1, but do_stars==0, so doing a multi-pass gas image
            ## 
            ##  -- many experiments here: doing gas isosurfaces with broad kernels
            ##       and overlaying a custom set of color tables after the fact seems
            ##       best. adjust opacity (kappa_units) and kernel width as needed. 
            ##       also, can set 'dynrange=0' to automatically produce a given 
            ##       fraction of saturated/unsaturated pixels
            ##
            gas_map_temperature_cuts=np.array([300., 2.0e4, 3.0e5 ])
            kernel_widths=np.array([0.8,0.3,0.6])
            #kernel_widths=np.array([0.7,0.3,0.7])

            # layer below is for zoom_dw large-scale GAS sims ONLY
            #gas_map_temperature_cuts=np.array([300., 1.0e4, 1.0e5 ])
            #kernel_widths=np.array([0.5,0.25,0.6])

            out_gas,out_u,out_g,out_r = rayproj.gas_raytrace_temperature( \
                gas_map_temperature_cuts, \
                x[ok], y[ok], z[ok], color_weights[ok], weights[ok], h_all[ok], \
                xrange=xr, yrange=yr, zrange=zr, pixels=pixels, \
                isosurfaces = 1, kernel_width=kernel_widths, \
                add_temperature_weights = 0 , KAPPA_UNITS = 2.0885*np.array([1.1,2.0,1.5]));
                #add_temperature_weights = 0 , KAPPA_UNITS = 2.0885*np.array([4.1,2.0,2.0]));
                
            image24, massmap, image_max, image_min = \
                makethreepic.make_threeband_image_process_bandmaps( out_r,out_g,out_u, \
                maxden=maxden,dynrange=dynrange,pixels=pixels, \
                color_scheme_nasa=nasa_colors,color_scheme_sdss=sdss_colors);

            image24 = makethreepic.layer_band_images(image24, massmap);                

    else:   ## threecolor==0    ( not a three-color image )
    
        if (do_with_colors==1):
            ##
            ## here we have a simple two-pass image where surface brightness is 
            ##   luminosity and a 'temperature' is color-encoded
            ##
            if (log_temp_wt==1): 
                color_weights=np.log10(color_weights); temp_min=np.log10(temp_min); temp_max=np.log10(temp_max);

            massmap_singlelayer,pic_singlelayer,massmap,image24 = \
            cmakepic.simple_makepic(x[ok],y[ok],weights=weights[ok],hsml=h_all[ok],\
                xrange=xr,yrange=yr,
                set_dynrng=dynrange,set_maxden=maxden,
                set_percent_maxden=set_percent_maxden,set_percent_minden=set_percent_minden, 
                color_temperature=1,temp_weights=color_weights[ok],
                pixels=pixels,invert_colorscale=invert_colors,
                set_temp_max=temp_max,set_temp_min=temp_min);

        else: ## single color-scale image
            ##
            ## and here, a single-color scale image for projected density of the quantity
            ##
            massmap,image_singledepth = \
            cmakepic.simple_makepic(x[ok],y[ok],weights=weights[ok],hsml=h_all[ok],\
                xrange=xr,yrange=yr,
                set_dynrng=dynrange,set_maxden=maxden,
                set_percent_maxden=set_percent_maxden,set_percent_minden=set_percent_minden, 
                color_temperature=0,pixels=pixels,invert_colorscale=1-invert_colors);

            ## convert to actual image using a colortable:
            my_cmap=matplotlib.cm.get_cmap('hot');
            image24 = my_cmap(image_singledepth/255.);

    ##
    ## whichever sub-routine you went to, you should now have a massmap (or set of them) 
    ##   and processed rgb 'image24' (which should always be re-makable from the massmaps)
    ##
    ## optionally can have great fun with some lighting effects to give depth:
    ##   (looks great, with some tuning)
    ##
    if (include_lighting==1):
        #light = matplotlib.colors.LightSource(azdeg=0,altdeg=65)
        light = viscolors.CustomLightSource(azdeg=0,altdeg=65)
        if (len(massmap.shape)>2):
            ## do some clipping to regulate the lighting:
            elevation = massmap.sum(axis=2)
            minden = maxden / dynrange
            elevation = (elevation - minden) / (maxden - minden)
            elevation[elevation < 0.] = 0.
            elevation[elevation > 1.] = 1.
            elevation *= maxden
            grad_max = maxden / 5.
            grad_max = maxden / 6.
            #image24_lit = light.shade_rgb(image24, massmap.sum(axis=2))
            image24_lit = light.shade_rgb(image24, elevation, vmin=-grad_max, vmax=grad_max)
        else:
            image24_lit = light.shade(image24, massmap)
        image24 = image24_lit
    
    plt.imshow(image24,origin='lower',interpolation='bicubic',aspect='normal');

    if ((add_gas_layer_to_image=='')==False):
        viscolors.load_my_custom_color_tables();
        plt.imshow(image_singledepth_extra_layer,origin='lower',interpolation='bicubic',aspect='normal',
            cmap=set_added_layer_ctable,alpha=set_added_layer_alpha)

    ## slap on the figure labels/scale bars
    if(show_scale_label==1): overlay_scale_label(xr_0,yr_0,ax_fig,c0='w')
    if(show_time_label==1): overlay_time_label(time,\
        ax_fig,cosmo=cosmo,c0='w',n_sig_add=0,tunit_suffix='Gyr')    
    
    filename=fname_base
    if(filename_set_manually!=''): filename=filename_set_manually
    plt.subplots_adjust(left=0.0,bottom=0.0,right=1-0.0,top=1-0.0,hspace=0.0,wspace=0.0);
    frame1.axes.xaxis.set_ticklabels([]) #no tick names
    frame1.axes.yaxis.set_ticklabels([]) #no tick names
    frame1.axes.get_xaxis().set_ticks([]) #no ticks
    frame1.axes.get_yaxis().set_ticks([]) #no ticks
    plt.tight_layout(pad=0,w_pad=0,h_pad=0)
    plt.savefig(filename+format,dpi=pixels,bbox_inches='tight',pad_inches=0)
    plt.close('all')

    ## save them:
    outfi = h5py.File(filename+'.dat','w')
    dset_im = outfi.create_dataset('image24',data=image24)
    dset_mm = outfi.create_dataset('massmap',data=massmap)
    dset_xr = outfi.create_dataset('xrange',data=xr_0)
    dset_yr = outfi.create_dataset('yrange',data=yr_0)
    dset_time = outfi.create_dataset('time',data=time)
    dset_cosmo = outfi.create_dataset('cosmological',data=cosmo)
    if ((add_gas_layer_to_image=='')==False):
        dset_im2 = outfi.create_dataset('image24_extralayer',data=image_singledepth_extra_layer)
    else:
        dset_im2 = outfi.create_dataset('image24_extralayer',data=np.array([-1]))
    outfi.close()
    gc.collect()
    return image24, massmap;


def xcompile(x_grid_previous,x_final,dlnx_previous,n):
    x_0 = x_grid_previous[-1]
    x_f = x_final
    dx_desired = np.log(x_f/x_0) / n
    lnx_grid_0 = np.log(x_grid_previous)
    lnx_new, dlnx_new = tcompile(lnx_grid_0,dx_desired,dlnx_previous,n)
    gc.collect()
    return np.exp(lnx_new), dlnx_new

def tcompile(t_grid,dt_desired,dt_previous,n):
    f_transition = 0.1
    dt_0 = dt_previous
    dt_1 = (dt_desired - 0.5*f_transition*dt_previous)/(1.0 - 0.5*f_transition)
    dt = np.zeros(np.round(n).astype(int))
    ng = 1.*np.arange(1,n+1)
    ok = (ng <= f_transition*n)
    dt[ok] = dt_0 + (dt_1-dt_0) * (1.*ng[ok])/(f_transition*n)
    ok = (ng > f_transition*n)
    dt[ok] = dt_1
    tg = np.cumsum(dt)
    t_grid_new = np.concatenate([t_grid, t_grid[-1] + tg])
    gc.collect()
    return t_grid_new, dt_1



    


    




def build_snapshot_list(snapdir , four_char=0):
    i=0; imax=750; snums=[-1];
    while (i < imax):
        fname,fname_base,fname_ext = check_if_filename_exists(snapdir,i,four_char=four_char)
        if(fname!='NULL'):
            snums=np.concatenate((snums,np.array([i])));
            if (i > imax-50): imax += 50;
        i += 1;
    print(snums)
    snums=snums[snums>=0];
    gc.collect()
    return snums






    

   
    

# builds the list of times to dump frames
def build_time_frame_grid( snapdir, snapshot_list, frames_per_gyr, \
        time_min=0, time_max=0, use_h0=1, four_char=0, cosmological=0, skip_bh=1 ):

    ## set times for frames
    if(time_min==0):
        PPP_head = gadget.readsnap( snapdir, snapshot_list[0], 1, header_only=1, \
            h0=use_h0,four_char=four_char,cosmological=cosmological,skip_bh=skip_bh )
        time_min = PPP_head['time'] ## set to first time in snapshot list
    if(time_max==0):
        PPP_head = gadget.readsnap( snapdir, snapshot_list[-1], 1, header_only=1, \
            h0=use_h0,four_char=four_char,cosmological=cosmological,skip_bh=skip_bh )
        time_max = PPP_head['time'] ## set to last time in snapshot list
    if(cosmological==1):
        time_min = cosmological_time(time_min) ## input in 'a_scale', not time (Gyr)
        time_max = cosmological_time(time_max)

    time_frame_grid = np.arange(time_min, time_max, 1./np.float(frames_per_gyr))

    if (cosmological==0): 
        a_scale_grid = time_frame_grid ## to avoid error below
    else:
        a_tmp = 10.**np.arange(-3.,0.001,0.001)
        t_tmp = cosmological_time(a_tmp)
        a_scale_grid = np.exp(np.interp(np.log(time_frame_grid),np.log(t_tmp),np.log(a_tmp)))
        a_scale_grid[a_scale_grid < a_tmp[0]] = a_tmp[0]

    gc.collect()
    return time_frame_grid, a_scale_grid



def checklen(x):
    gc.collect()
    return len(np.array(x,ndmin=1));


def ok_scan(input,xmax=1.0e10,pos=0):
    if (pos==1):
        gc.collect()
        return (np.isnan(input)==False) & (np.abs(input)<=xmax) & (input > 0.);
    else:
        gc.collect()
        return (np.isnan(input)==False) & (np.abs(input)<=xmax);


def snap_ext(snum,four_char=0):
	ext='00'+str(snum);
	if (snum>=10): ext='0'+str(snum)
	if (snum>=100): ext=str(snum)
	if (four_char==1): ext='0'+ext
	if (snum>=1000): ext=str(snum)
	gc.collect()
	return ext;

	
def get_gas_xray_luminosity(ppp):
    brems = gadget.gas_xray_brems( \
        ppp['m'], ppp['u'], ppp['rho'], ppp['ne'], ppp['nh'] );
    gc.collect()
    return brems; ## can also add metal line cooling luminosities later

    

def overlay_scale_label(xr,yr,figure_axis,c0='w'):
    ddx=0.5*(xr[1]-xr[0]);
    if (ddx>0.00002): dx=0.00001; dxl='0.01 pc';
    if (ddx>0.0002): dx=0.0001; dxl='0.1 pc';
    if (ddx>0.002): dx=0.001; dxl='1 pc';
    if (ddx>0.02): dx=0.01; dxl='10 pc';
    if (ddx>0.2): dx=0.1; dxl='100 pc';
    if (ddx>2.): dx=1.; dxl='1 kpc';
    if (ddx>20.): dx=10.; dxl='10 kpc';
    if (ddx>200.): dx=100.; dxl='100 kpc';
    if (ddx>2000.): dx=1000.; dxl='1 Mpc';
    if (ddx>20000.): dx=10000.; dxl='10 Mpc';
    
    xlen=(xr[1]-xr[0])
    ylen=(yr[1]-yr[0])
    xoff = (0.25+0.02)*ddx / xlen
    yoff = 0.025*(yr[1]-yr[0]) / ylen
    xr_new = np.array([xoff-dx/xlen*0.5,xoff+dx/xlen*0.5])
    yr_new = np.array([yoff,yoff])
    
    plt.text(xoff,1.75*yoff,dxl,color=c0,\
        horizontalalignment='center',verticalalignment='baseline',\
        transform=figure_axis.transAxes,fontsize=3)
    figure_axis.autoscale(enable=False,axis='both',tight=True)
    figure_axis.plot(xr_new,yr_new,color=c0,linewidth=0.7,transform=figure_axis.transAxes)

    
def overlay_time_label(time, figure_axis, cosmo=0, c0='w', n_sig_add=0, tunit_suffix='Gyr'):
    if(cosmo==0): 
        time_to_use = time ## absolute time in Gyr
        prefix = ''
        suffix = tunit
    else: ## 
        time_to_use = 1./time - 1. ## redshift
        prefix = 'z='
        suffix = ''
        
    n_sig = 2
    if time_to_use >= 10: n_sig += 0
    if time_to_use < 1.: n_sig += 1
    t_str = round_to_n( time_to_use, n_sig+n_sig_add )
    label_str = prefix+t_str+suffix
    xoff=0.03; yoff=1.-0.025;
    plt.text(xoff,yoff,label_str,color=c0,\
        horizontalalignment='left',verticalalignment='top',\
        transform=figure_axis.transAxes,fontsize=3)


def round_to_n(x, n):
    ''' Utility function used to round labels to significant figures for display purposes  from: http://mail.python.org/pipermail/tutor/2004-July/030324.html'''
    if n < 1:
        raise ValueError("number of significant digits must be >= 1")

    # show everything as floats (preference; can switch using code below to showing eN instead
    format = "%." +str(n-1) +"f"
    as_string=format % x
    gc.collect()
    return as_string

    # Use %e format to get the n most significant digits, as a string.
    format = "%." + str(n-1) + "e"
    as_string = format % x
    if as_string[-3:] in ['+00', '+01', '+02', '+03','-01', '-02', '-03']:
        #then number is 'small', show this as a float
        format = "%." +str(n-1) +"f"
        as_string=format % x
    gc.collect()
    return as_string




##############################


## procedure which will return, for a given input vector A_in, 
##    the perpendicular unit vectors B_out and C_out 
##    which form perpendicular axes to A
def return_perp_vectors(a, LOUD=0):
    eps = 1.0e-10
    a = np.array(a,dtype='f');
    a /= np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
    for i in range(len(a)):
        if (a[i]==0.): a[i]=eps;
        if (a[i]>=1.): a[i]=1.-eps;
        if (a[i]<=-1.): a[i]=-1.+eps;
    a /= np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
    ax=a[0]; ay=a[1]; az=a[2];

    ## use a fixed rotation of the a-vector by 90 degrees:
    ## (this anchors the solution so it changes *continously* as a changes)
    t0=np.double(np.pi/2.e0);
    bx=0.*ax; by=np.cos(t0)*ay-np.sin(t0)*az; bz=np.sin(t0)*ay+np.cos(t0)*az;
    ## c-sign is degenerate even for 'well-chosen' a and b: gaurantee right-hand 
    ##   rule is obeyed by defining c as the cross product: a x b = c
    cx=(ay*bz-az*by); cy=-(ax*bz-az*bx); cz=(ax*by-ay*bx); 
    B_out=np.zeros(3); C_out=np.zeros(3);
    B_out[:]=[bx,by,bz]; C_out[:]=[cx,cy,cz];
    
    if (LOUD==1):
        print(a)
        print(B_out)
        print(C_out)
        print('a_tot=',ax*ax+ay*ay+az*az)
        print('b_tot=',bx*bx+by*by+bz*bz)
        print('c_tot=',cx*cx+cy*cy+cz*cz)
        print('adotb=',ax*bx+ay*by+az*bz)
        print('adotc=',ax*cx+ay*cy+az*cz)
        print('bdotc=',bx*cx+by*cy+bz*cz)
    return B_out, C_out


def frame_ext(snum,four_char=1):
	ext='00'+str(snum);
	if (snum>=10): ext='0'+str(snum)
	if (snum>=100): ext=str(snum)
	if (four_char==1): ext='0'+ext
	if (snum>=1000): ext=str(snum)
	gc.collect()
	return ext;



def ok_scan(input,xmax=1.0e10,pos=0):
    if (pos==1):
        gc.collect()
        return (np.isnan(input)==False) & (np.abs(input)<=xmax) & (input > 0.);
    else:
        gc.collect()
        return (np.isnan(input)==False) & (np.abs(input)<=xmax);


def idmatch(id_i, id_f): 
    ## match items in set (i) to those in set (f)
    ##  --- this is particularly speedy when the ids are all unique (case here) --
    index_dict_i = dict((k,i) for i,k in enumerate(id_i));
    index_dict_f = dict((k,i) for i,k in enumerate(id_f));
    inter = set(id_i).intersection(set(id_f));
    indices_i = np.array([ index_dict_i[x] for x in inter ]);
    indices_f = np.array([ index_dict_f[x] for x in inter ]);
    gc.collect()
    return indices_i, indices_f;


def compile_matched_ids( id_i, id_f ):
    sub_id_i, sub_id_f = idmatch(id_i,id_f)
    
    nomatch_from_i_in_f = (id_i > -1.0e40) ## should return all true
    if (sub_id_i.size > 0):
        nomatch_from_i_in_f[sub_id_i] = False ## is matched
    
    nomatch_from_f_in_i = (id_f > -1.0e40) 
    if (sub_id_f.size > 0):
        nomatch_from_f_in_i[sub_id_f] = False
    
    gc.collect()
    return sub_id_i, sub_id_f, nomatch_from_i_in_f, nomatch_from_f_in_i

    
def cross(x,y):
    #return np.cross(x,y,axis=1)
    c=0.*x
    c[:,0] = x[:,1]*y[:,2] - x[:,2]*y[:,1]
    c[:,1] =-x[:,0]*y[:,2] + x[:,2]*y[:,0]
    c[:,2] = x[:,0]*y[:,1] - x[:,1]*y[:,0]
    gc.collect()
    return c


def cosmological_time(a,h=0.71,Omega_M=0.27):
    ## exact solution for a flat universe
    x=Omega_M/(1.-Omega_M) / (a*a*a);
    t=(2./(3.*np.sqrt(1.-Omega_M))) * np.log( np.sqrt(x) / (-1. + np.sqrt(1.+x)) );
    t *= 13.777 * (0.71/h); ## in Gyr
    gc.collect()
    return t;
    
    
def fcor(x):
    gc.collect()
    return np.array(x,dtype='f',ndmin=1)

def vfloat(x):
    gc.collect()
    return x.ctypes.data_as(ctypes.POINTER(ctypes.c_float));


##############################################


def interpolation_for_movie(dt, delta_t, \
	    m_i,h_i,c_i,zm_i,\
	    m_f,h_f,c_f,zm_f, \
	    x_in,y_in,z_in,vx_in,vy_in,vz_in, \
	    x_fn,y_fn,z_fn,vx_fn,vy_fn,vz_fn, \
	    s_i,s_f,nm_i_f,nm_f_i, \
	    cen_i,cen_f , use_polar_interpolation=1, 
	    polar_vs_cylindrical = True, 
	    set_dx_frame = 0., stars = 0):
    
    x_i = np.copy(x_in) - cen_i[0]
    y_i = np.copy(y_in) - cen_i[1]
    z_i = np.copy(z_in) - cen_i[2]
    vx_i = np.copy(vx_in) - (cen_f[0]-cen_i[0])/delta_t
    vy_i = np.copy(vy_in) - (cen_f[1]-cen_i[1])/delta_t
    vz_i = np.copy(vz_in) - (cen_f[2]-cen_i[2])/delta_t
    x_f = np.copy(x_fn) - cen_f[0]
    y_f = np.copy(y_fn) - cen_f[1]
    z_f = np.copy(z_fn) - cen_f[2]
    vx_f = np.copy(vx_fn) - (cen_f[0]-cen_i[0])/delta_t
    vy_f = np.copy(vy_fn) - (cen_f[1]-cen_i[1])/delta_t
    vz_f = np.copy(vz_fn) - (cen_f[2]-cen_i[2])/delta_t
    
    #for ui,uin,j in zip([x_i,y_i,z_i],[x_in,y_in,z_in],[0,1,2]): ui = uin - cen_i[j]
    #for uf,ufn,j in zip([x_f,y_f,z_f],[x_fn,y_fn,z_fn],[0,1,2]): uf = ufn - cen_f[j]
    #for vui,vuin,j in zip([vx_i,vy_i,vz_i],[vx_in,vy_in,vz_in],[0,1,2]): vui = vuin - (cen_f[j]-cen_i[j])/delta_t
    #for vuf,vufn,j in zip([vx_f,vy_f,vz_f],[vx_fn,vy_fn,vz_fn],[0,1,2]): vuf = vufn - (cen_f[j]-cen_i[j])/delta_t
    pos_interp_flag = np.array([0.])
    if(stars!=0): 
        r_i = np.sqrt(x_i[s_i]*x_i[s_i] + y_i[s_i]*y_i[s_i] + z_i[s_i]*z_i[s_i])
        r_f = np.sqrt(x_f[s_f]*x_f[s_f] + y_f[s_f]*y_f[s_f] + z_f[s_f]*z_f[s_f])
        d_r = np.sqrt((x_i[s_i]-x_f[s_f])**2 + (y_i[s_i]-y_f[s_f])**2 + (z_i[s_i]-z_f[s_f])**2)
        r_eff = 1. / (1./r_i + 1./r_f)
        dvx=(x_f[s_f]-x_i[s_i])/delta_t; dvy=(y_f[s_f]-y_i[s_i])/delta_t; dvz=(z_f[s_f]-z_i[s_i])/delta_t
        dv2_i=(vx_i[s_i]-dvx)**2 + (vy_i[s_i]-dvy)**2 + (vz_i[s_i]-dvz)**2
        dv2_f=(vx_f[s_f]-dvx)**2 + (vy_f[s_f]-dvy)**2 + (vz_f[s_f]-dvz)**2
        dv = np.sqrt(dv2_i + dv2_f); v0 = np.sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
        pos_interp_flag = (dv/v0 > 1.0) & (dv > 100.) & (d_r/r_eff > 2.2) & (d_r > 0.01*set_dx_frame) & (r_eff < 0.2*set_dx_frame)

    if (use_polar_interpolation==1):
        # convert to cylindrical coordinates, with the axis defined as cumulative ang-mom axis inside each orbit
        r_i,vr_i,zc_i,vzc_i,phi_i,vphi_i,j_hat_i = xyz_to_cylindrical(x_i,y_i,z_i,vx_i,vy_i,vz_i)
        r_f,vr_f,zc_f,vzc_f,phi_f,vphi_f,j_hat_f = xyz_to_cylindrical(x_f,y_f,z_f,vx_f,vy_f,vz_f)
        if(stars!=0): 
            z0 = np.sqrt(zc_i[s_i]*zc_i[s_i] + zc_f[s_f]*zc_f[s_f])
            dz = np.abs(zc_i[s_i]-zc_f[s_f])
            #pos_interp_flag = pos_interp_flag & (z0/r_eff > 0.5) & (z0 > 0.5)
            #pos_interp_flag = pos_interp_flag & (dz/z0 > 1.0) & (r_eff > 0.5) & (r_eff < 3.) 
            pos_interp_flag = pos_interp_flag & (r_eff < 3.) & (dz > 0.1) & (dz/z0 > 1.) & (z0 > 0.1)
            
            
            #pos_interp_flag = pos_interp_flag & (z0/r_eff > 0.5) & (z0 > 0.25) & (dz/z0 > 0.75) & (r_eff < 3.)
        if(polar_vs_cylindrical == True):
            # convert to polar coordinates, where z points out of plane of orbit for each particle
            r_i,vr_i,phi_i,vphi_i,j_hat_i = xyz_to_polar(x_i,y_i,z_i,vx_i,vy_i,vz_i)
            r_f,vr_f,phi_f,vphi_f,j_hat_f = xyz_to_polar(x_f,y_f,z_f,vx_f,vy_f,vz_f)
        
        # r gets treated like a normal spatial coordinate (except cannot go <0; ignore for now?)
        r = interpolation_for_movie_pos(dt,delta_t, r_i,r_f,vr_i,vr_f, s_i,s_f,nm_i_f,nm_f_i,
            order=1,pos_interp_flag=pos_interp_flag,threshold_for_overshoot=1.e-10)
        # predict the 'final' phi, theta, knowing that these are actually periodic:
        # note for phi, need third-order: inexact matching can leave 'wedges' missing from orbits
        phi = interpolation_for_movie_pos(dt,delta_t, phi_i,phi_f,
            vphi_i,vphi_f, s_i,s_f,nm_i_f,nm_f_i, periodic=1,order=3,
            threshold_for_overshoot=0.25)
        x_t = r * np.cos(phi);  y_t = r * np.sin(phi);
        j_hat = interpolation_for_movie_jhat(dt,delta_t, j_hat_i,j_hat_f, s_i,s_f,nm_i_f,nm_f_i )

        #if(stars!=0): 
        #    pos_interp_flag = (r < 

        ## define the absolute 'phi' relative to an arbitrary (but fixed) 90-degree rotation of the angular momentum vector j
        x_jhat = 0.*j_hat
        x_jhat[:,0]=0.*j_hat[:,0]; x_jhat[:,1]=-j_hat[:,2]; x_jhat[:,2]=j_hat[:,1]
        jj = np.sqrt( x_jhat[:,0]*x_jhat[:,0] + x_jhat[:,1]*x_jhat[:,1] + x_jhat[:,2]*x_jhat[:,2] )
        for j in [0,1,2]: x_jhat[:,j] /= jj
        ## generate y-vector by cross-product a x b = c (gaurantees right-hand rule)
        y_jhat = cross( j_hat, x_jhat )
        jj = np.sqrt( y_jhat[:,0]*y_jhat[:,0] + y_jhat[:,1]*y_jhat[:,1] + y_jhat[:,2]*y_jhat[:,2] )
        for j in [0,1,2]: y_jhat[:,j] /= jj
        z_t = 0.0 * x_t
        if(polar_vs_cylindrical == False):
            z_t = interpolation_for_movie_pos(dt,delta_t, zc_i,zc_f,vzc_i,vzc_f, s_i,s_f,nm_i_f,nm_f_i,pos_interp_flag=pos_interp_flag)#,order=1)
        ## now project back into our lab coordinate system
        x = x_t*x_jhat[:,0] + y_t*y_jhat[:,0] + z_t*j_hat[:,0]
        y = x_t*x_jhat[:,1] + y_t*y_jhat[:,1] + z_t*j_hat[:,1]
        z = x_t*x_jhat[:,2] + y_t*y_jhat[:,2] + z_t*j_hat[:,2]

        ## use cartesian interpolation for some points
        x_c=interpolation_for_movie_pos(dt,delta_t, x_i,x_f,vx_i,vx_f, s_i,s_f,nm_i_f,nm_f_i,pos_interp_flag=pos_interp_flag)#,order=1)
        y_c=interpolation_for_movie_pos(dt,delta_t, y_i,y_f,vy_i,vy_f, s_i,s_f,nm_i_f,nm_f_i,pos_interp_flag=pos_interp_flag)#,order=1)
        z_c=interpolation_for_movie_pos(dt,delta_t, z_i,z_f,vz_i,vz_f, s_i,s_f,nm_i_f,nm_f_i,pos_interp_flag=pos_interp_flag)#,order=1)

        ## combine velocities into initial and final matched vectors
        Nmatched=x_i[s_i].size;
        Nunmatched_i=x_i[nm_i_f].size;
        Nunmatched_f=x_f[nm_f_i].size;
        u=np.zeros(Nmatched+Nunmatched_i+Nunmatched_f);
        z0=0.*u+1.e-10; dz=0.*u; vr=0.*u; vphi=0.*u; rr=0.*u; dr=0.*u; jj=0.*u; vz=0.*u; rvphi_i=r_i*vphi_i; rvphi_f=r_f*vphi_f;
        if (Nmatched>0):
            vr[0:Nmatched] = np.sqrt((vr_i[s_i]**2.+vr_f[s_f]**2.)/2.)
            vphi[0:Nmatched] = np.sqrt((rvphi_i[s_i]**2.+rvphi_f[s_f]**2.)/2.)
            rr[0:Nmatched] = np.sqrt((r_i[s_i]**2.+r_f[s_f]**2.)/2.)
            dr[0:Nmatched] = np.abs(r_i[s_i]-r_f[s_f]) / np.abs(r_i[s_i]+r_f[s_f])
            jj[0:Nmatched] = (j_hat_i[s_i,0]*j_hat_f[s_f,0] + j_hat_i[s_i,1]*j_hat_f[s_f,1] + j_hat_i[s_i,2]*j_hat_f[s_f,2]) / \
                np.sqrt((j_hat_i[s_i,0]*j_hat_i[s_i,0] + j_hat_i[s_i,1]*j_hat_i[s_i,1] + j_hat_i[s_i,2]*j_hat_i[s_i,2]) * \
                        (j_hat_f[s_f,0]*j_hat_f[s_f,0] + j_hat_f[s_f,1]*j_hat_f[s_f,1] + j_hat_f[s_f,2]*j_hat_f[s_f,2]))
            vz[0:Nmatched] = 0
            if(use_polar_interpolation==1):
                z0[0:Nmatched] = np.sqrt(zc_i[s_i]*zc_i[s_i] + zc_f[s_f]*zc_f[s_f])
                dz[0:Nmatched] = np.abs(zc_i[s_i] - zc_f[s_f])
            if(polar_vs_cylindrical == False):
                vz[0:Nmatched] = np.sqrt((vzc_i[s_i]**2.+vzc_f[s_f]**2.)/2.)
        if(Nunmatched_i>0):
            vr[Nmatched:Nmatched+Nunmatched_i] = vr_i[nm_i_f]
            vphi[Nmatched:Nmatched+Nunmatched_i] = rvphi_i[nm_i_f]
            rr[Nmatched:Nmatched+Nunmatched_i] = r_i[nm_i_f]
            rr[Nmatched:Nmatched+Nunmatched_i] = -1
            jj[Nmatched:Nmatched+Nunmatched_i] = 1
            vz[Nmatched:Nmatched+Nunmatched_i] = 1
        if(Nunmatched_f>0):
            vr[Nmatched+Nunmatched_i:Nmatched+Nunmatched_i+Nunmatched_f] = vr_f[nm_f_i]
            vphi[Nmatched+Nunmatched_i:Nmatched+Nunmatched_i+Nunmatched_f] = rvphi_f[nm_f_i]
            rr[Nmatched+Nunmatched_i:Nmatched+Nunmatched_i+Nunmatched_f] = r_f[nm_f_i]
            rr[Nmatched+Nunmatched_i:Nmatched+Nunmatched_i+Nunmatched_f] = -1
            jj[Nmatched+Nunmatched_i:Nmatched+Nunmatched_i+Nunmatched_f] = 1
            vz[Nmatched+Nunmatched_i:Nmatched+Nunmatched_i+Nunmatched_f] = 1
        ## attempt to guess what's 'rotation supported' and not, and use appropriate interpolation
        vr2 = vr*vr + vz*vz
        v2 = vr2 + vphi*vphi
        
        ## good on large scales
        no_rot_support = (z0>2.) | ( z0/rr > 0.6 ) | ( (dz/z0>2.)&(z0 > 1.) ) | ( rr == -1 ) | (jj < 0.8) | ( vr2 > v2/2. ) | ( rr > 12. ) # usually does ok, too much 'breathing' at late times
        no_rot_support_favorcartesian = no_rot_support
        
        ## good on small scales
        no_rot_support = (rr == -1) | (rr > 12.) | (set_dx_frame>60.) | ((dr/rr > 2.2)&(rr>2)&(set_dx_frame>30.)) ##| (jj < (-1.+2.*set_dx_frame/60.))
        favorcartesian = (np.random.RandomState(m_i.size+m_f.size).rand(rr.size) < np.exp(-(30./set_dx_frame)**2.))
        no_rot_support[favorcartesian] = no_rot_support_favorcartesian[favorcartesian]
        
        x[no_rot_support] = x_c[no_rot_support]
        y[no_rot_support] = y_c[no_rot_support]
        z[no_rot_support] = z_c[no_rot_support]
    else:
        x=interpolation_for_movie_pos(dt,delta_t, x_i,x_f,vx_i,vx_f, s_i,s_f,nm_i_f,nm_f_i,pos_interp_flag=pos_interp_flag)
        y=interpolation_for_movie_pos(dt,delta_t, y_i,y_f,vy_i,vy_f, s_i,s_f,nm_i_f,nm_f_i,pos_interp_flag=pos_interp_flag)
        z=interpolation_for_movie_pos(dt,delta_t, z_i,z_f,vz_i,vz_f, s_i,s_f,nm_i_f,nm_f_i,pos_interp_flag=pos_interp_flag)
                

    x = x + (cen_f[0]-cen_i[0]) * dt / delta_t + cen_i[0]
    y = y + (cen_f[1]-cen_i[1]) * dt / delta_t + cen_i[1]
    z = z + (cen_f[2]-cen_i[2]) * dt / delta_t + cen_i[2]

    h=interpolation_for_movie_scalar(dt,delta_t, h_i,h_f, s_i,s_f,nm_i_f,nm_f_i, hsml_rescale=1, pos_interp_flag=pos_interp_flag)
    m=interpolation_for_movie_scalar(dt,delta_t, m_i,m_f, s_i,s_f,nm_i_f,nm_f_i, mass_rescale=1, pos_interp_flag=pos_interp_flag)
    c=interpolation_for_movie_scalar(dt,delta_t, c_i,c_f, s_i,s_f,nm_i_f,nm_f_i)
    #c=np.log10(interpolation_for_movie_scalar(dt,delta_t, 10.**c_i,10.**c_f, s_i,s_f,nm_i_f,nm_f_i))
    zm=interpolation_for_movie_scalar(dt,delta_t, zm_i,zm_f, s_i,s_f,nm_i_f,nm_f_i)

    gc.collect()
    return x,y,z,m,h,c,zm;






def interpolation_for_movie_pos( dt, delta_t, u_i,u_f, vu_i,vu_f, \
    sub_id_i,sub_id_f, nomatch_from_i_in_f,nomatch_from_f_in_i , periodic=0, order=3,
    pos_interp_flag=np.array([0.]), threshold_for_overshoot=0.25):

    #order = 1

    dt=np.array(dt); delta_t=np.array(delta_t); ## make sure casting is ok
    tau = dt/delta_t;
    Nmatched=u_i[sub_id_i].size;
    Nunmatched_i=u_i[nomatch_from_i_in_f].size;
    Nunmatched_f=u_f[nomatch_from_f_in_i].size;
    u=np.zeros(Nmatched+Nunmatched_i+Nunmatched_f);

    # first compute for objects in (i) with a match in (f)
    if (Nmatched>0):
        xi = np.copy(u_i[sub_id_i])  ; xf = np.copy(u_f[sub_id_f])  ;
        vi = np.copy(vu_i[sub_id_i]) ; vf = np.copy(vu_f[sub_id_f]) ;
        
        order_vec = 0.*xi + order
        order_tol = threshold_for_overshoot
        vdenom = 1. / (order_tol * (1.e-36 + np.abs(vf+vi)))
        xdenom = 1. / (order_tol * (1.e-36 + np.abs(xf+xi)))
        
        order_vec[(np.abs(vf-vi)*vdenom  > 1)] = 1
        
        if (periodic==1):
            ## u is periodic in 2pi, so estimate its final position and wrap coordinate
            ## guess where its going to be to get 'best estimate' of the number of revolutions
            #x_exp = xi + 0.5*(vf+vi)*delta_t # actually can get some errors when vf/i switch signs...
            x_exp = xi + vi*delta_t
            ## now pin the final location to the appropriate value
            xf = get_closest_periodic_value( x_exp, xf , xi)
        
        ## check for over/under-shooting, and if present, lower the reconstruction order
        #if(order==3):
        x0=xi; x1=vi*delta_t; x2=3.*(xf-xi)-(2.*vi+vf)*delta_t; x3=-2.*(xf-xi)+(vi+vf)*delta_t  
        # velocity extremum 
        t_x = -x2/(3.*x3 + 1.e-40)
        v_x=(x1-x2*t_x)/delta_t; 
        dv = np.minimum(np.abs(v_x-vi) , np.abs(v_x-vf))
        ok_tmp = (dv*vdenom > 1) & (t_x > 0.) & (t_x < 1.)
        order_vec[ok_tmp] = np.minimum(order_vec[ok_tmp],2 + np.zeros(order_vec[ok_tmp].size))
        # position extremum 1 (b/c quadratic solution)
        q=x2*x2 - 3.*x1*x3
        qq=0.*q; 
        ok_t=(q>0); qq[ok_t] = t_x[ok_t] + np.sqrt(q[ok_t])/(3.*x3[ok_t])
        ok = (q > 0) & (qq > 0.) & (qq < 1.)
        xx=x0 + x1*qq + x2*qq*qq + x3*qq*qq*qq
        dx = np.minimum(np.abs(xx-xi) , np.abs(xx-xf))
        ok_tmp = (dx*xdenom > 1) & ok
        ok_tmp = (q > 0) & (qq > 0.) & (qq < 1.)
        order_vec[ok_tmp] = np.minimum(order_vec[ok_tmp],2 + np.zeros(order_vec[ok_tmp].size))
        # position extremum 2 (b/c quadratic solution)
        qq=0.*q; 
        ok_t=(q>0); qq[ok_t] = t_x[ok_t] - np.sqrt(q[ok_t])/(3.*x3[ok_t])
        ok = (q > 0) & (qq > 0.) & (qq < 1.)
        xx=x0 + x1*qq + x2*qq*qq + x3*qq*qq*qq
        dx = np.minimum(np.abs(xx-xi) , np.abs(xx-xf))
        ok_tmp = (dx*xdenom > 1) & ok
        ok_tmp = (q > 0) & (qq > 0.) & (qq < 1.)
        order_vec[ok_tmp] = np.minimum(order_vec[ok_tmp],2 + np.zeros(order_vec[ok_tmp].size))

        #if(order==2):
        x0=xi; x2=(vf-vi)*delta_t/2.; x1=(xf-xi)-x2;
        t_x=-x1/(2.*x2 + 1.e-40)
        xx = x0 + (x1/2.)*t_x;
        dx = np.minimum(np.abs(xx-xi) , np.abs(xx-xf))
        ok_tmp = (dx*xdenom > 1) & (t_x > 0.) & (t_x < 1.) & (order_vec <= 2.)
        ok_tmp = (t_x > 0.) & (t_x < 1.) & (order_vec <= 2.)
        order_vec[ok_tmp] = np.minimum(order_vec[ok_tmp],1 + np.zeros(order_vec[ok_tmp].size))
        
        #order_vec = 0.*order_vec + 3.
        
        ## third-order interpolation: enables exact matching of x, v, but can 'overshoot'
        x0=xi; x1=vi*delta_t; x2=3.*(xf-xi)-(2.*vi+vf)*delta_t; x3=-2.*(xf-xi)+(vi+vf)*delta_t  
        u_o3 = x0 + x1*tau + x2*tau*tau + x3*tau*tau*tau
        ## second-order: minimizes absolute velocity difference (more stable, no exact velocity matching)
        x0=xi; x2=(vf-vi)*delta_t/2.; x1=(xf-xi)-x2;
        u_o2 = x0 + x1*tau + x2*tau*tau
        ## use linear interpolation below if for some reason this is unstable: least accurate, but most stable
        x0=xi; x1=(xf-xi)
        u_o1 = x0 + x1*tau
        
        #x0=xi; x1=vi*delta_t; x4=9.*(vf-vi)*delta_t/4.; 
        #x4=7.*(vf-vi)*delta_t/4.; 
        #x3=2.*(xi-(xf+x4))+(vf+vi)*delta_t; x2=0.5*((vf-vi)*delta_t-3.*x3-4.*x4);
        #u_o1 = x0 + x1*tau + x2*tau*tau + x3*tau*tau*tau + x4*tau*tau*tau*tau 
        
        #u_o1 = xi + (xf-xi)*tau + 0.5*((vf+vi)*delta_t+2.*(xi-xf))*np.sin(2.*np.pi*tau)/(2.*np.pi)
        
        u_tmp = u_o1
        order_vec=np.round(order_vec)
        ok = (order_vec == 3)
        u_tmp[ok] = u_o3[ok]
        ok = (order_vec == 2)
        u_tmp[ok] = u_o2[ok]
        #u[0:Nmatched] = u_tmp[0:Nmatched]
        
        #v0 = (xf-xi)/delta_t
        #dv = np.sqrt((vi-v0)**2 + (vf-v0)**2)
        #ok_tmp = (dv/(1.e-40 + np.abs(v0)) > 0.5) # too many at x_i, but does resolve central 'ring'
        #ok_tmp = (dv/(1.e-40 + np.abs(v0)) > 5.) # ring re-appears, potentially better?
        #ok_tmp = (dv/(1.e-40 + np.abs(v0)) > 1.6) # somewhere in-between, but not great
        #u_o1[ok_tmp] = xf[ok_tmp]
        #u_o1 = xf
        
        #dx = np.abs(xf-xi)
        #ok_tmp = (dx > 1.) # ok, good, but noisy
        #ok_tmp = (dx > 10.) # ring remains
        #ok_tmp = (dx > 3.) # ring remains
        
        if(pos_interp_flag.size > 1):
            ## smoothed step-function interpolation between snapshots
            width_ftau = 0.05
            width_ftau = 0.25
            q = 2.*np.pi/width_ftau
            f0 = 0.5-np.arctan(q*(0.-0.5))/np.pi
            f1 = 0.5-np.arctan(q*(1.-0.5))/np.pi
            f_tau = (0.5-np.arctan(q*(tau-0.5))/np.pi - f1) / (f0 - f1)
            sharp_interp = (xi + vi*dt) * f_tau + (xf + vf*(dt-delta_t)) * (1.-f_tau)
            #sharp_interp = (xi + 0*vi*dt) * f_tau + (xf + 0*vf*(dt-delta_t)) * (1.-f_tau)
            #u_tmp[pos_interp_flag] = sharp_interp[pos_interp_flag]
            #u_tmp[pos_interp_flag] = u_o3[pos_interp_flag]
            #if(tau < 0.5):
            #    u_tmp[pos_interp_flag] = xi[pos_interp_flag]
            #else:
            #    u_tmp[pos_interp_flag] = xf[pos_interp_flag]
            
        u[0:Nmatched] = u_tmp[0:Nmatched]
        
        
    ## now unmatched: first those in (i) with no match in (f)
    if(Nunmatched_i>0):
        u[Nmatched:Nmatched+Nunmatched_i] = \
            u_i[nomatch_from_i_in_f]+vu_i[nomatch_from_i_in_f]*dt;
	## now unmatched: now those in (f) with no match in (i)
    if(Nunmatched_f>0):
        u[Nmatched+Nunmatched_i:Nmatched+Nunmatched_i+Nunmatched_f] = \
            u_f[nomatch_from_f_in_i]+vu_f[nomatch_from_f_in_i]*(dt-delta_t);

    gc.collect()
    return u;


def interpolation_for_movie_jhat( dt, delta_t, j_i,j_f, \
    sub_id_i,sub_id_f, nomatch_from_i_in_f,nomatch_from_f_in_i ):

    dt=np.array(dt); delta_t=np.array(delta_t); ## make sure casting is ok
    tau = dt/delta_t;
    u_i = np.zeros(j_i[:,0].size)
    u_f = np.zeros(j_f[:,0].size)
    Nmatched=u_i[sub_id_i].size;
    Nunmatched_i=u_i[nomatch_from_i_in_f].size;
    Nunmatched_f=u_f[nomatch_from_f_in_i].size;
    j_hat=np.zeros((Nmatched+Nunmatched_i+Nunmatched_f,3));

    # first compute for objects in (i) with a match in (f)
    if (Nmatched>0):
        ## initialize matrices
        pt_i = np.zeros(Nmatched); pt_f = 0.*pt_i; 
        p_i = np.zeros((Nmatched,3)); p_f = 0.*p_i; v_f = 0.*p_f
        for j in [0,1,2]:
            p_i[:,j] = j_i[sub_id_i,j]
            p_f[:,j] = j_f[sub_id_f,j]
            pt_i += p_i[:,j]*p_i[:,j]
            pt_f += p_f[:,j]*p_f[:,j]
        ## normalize vectors to unity
        pt_i = np.sqrt(pt_i); pt_f = np.sqrt(pt_f);
        for j in [0,1,2]:
            p_i[:,j] /= pt_i
            p_f[:,j] /= pt_f
            
        ## build a rotation matrix between the j_hat values at two different timesteps:
        cos_rot_ang = p_i[:,0]*p_f[:,0] + p_i[:,1]*p_f[:,1] + p_i[:,2]*p_f[:,2]
        angle = np.arccos( cos_rot_ang )

        angle_to_rotate = angle * tau 
        cos_ang = np.cos(angle_to_rotate)
        sin_ang = np.sin(angle_to_rotate)

        ## build perpendicular vectors and normalize
        k_rot = cross( p_i, p_f )
        kk=0.*pt_i;
        for j in [0,1,2]: kk += k_rot[:,j]*k_rot[:,j]
        kk = np.sqrt(kk)
        for j in [0,1,2]: k_rot[:,j] /= kk
        
        k_cross_p = cross( k_rot, p_i )
        kk=0.*pt_i;
        for j in [0,1,2]: kk += k_cross_p[:,j]*k_cross_p[:,j]
        kk = np.sqrt(kk)
        for j in [0,1,2]: k_cross_p[:,j] /= kk
        
        k_dot_p = k_rot[:,0]*p_i[:,0] + k_rot[:,1]*p_i[:,1] + k_rot[:,2]*p_i[:,2]
        for j in [0,1,2]:
            v_f[:,j] = p_i[:,j]*cos_ang + k_cross_p[:,j]*sin_ang + k_rot[:,j]*k_dot_p*(1.-cos_ang)

        for j in [0,1,2]:
            j_hat[0:Nmatched,j] = v_f[:,j]  ## now have the j_hat vector for this time

    ## now unmatched: first those in (i) with no match in (f)
    if(Nunmatched_i>0):
        for j in [0,1,2]:
            j_hat[Nmatched:Nmatched+Nunmatched_i,j] = j_i[nomatch_from_i_in_f,j]

	## now unmatched: now those in (f) with no match in (i)
    if(Nunmatched_f>0):
        for j in [0,1,2]:
            j_hat[Nmatched+Nunmatched_i:Nmatched+Nunmatched_i+Nunmatched_f] = j_f[nomatch_from_f_in_i]

    gc.collect()
    return j_hat;



def interpolation_for_movie_scalar(dt, delta_t, u_i,u_f, \
    sub_id_i,sub_id_f, nomatch_from_i_in_f,nomatch_from_f_in_i, \
    mass_rescale=0, hsml_rescale=0, pos_interp_flag=np.array([0.]) ):

    dt=np.array(dt); delta_t=np.array(delta_t); ## make sure casting is ok
    tau = dt/delta_t;
    Nmatched=u_i[sub_id_i].size;
    Nunmatched_i=u_i[nomatch_from_i_in_f].size;
    Nunmatched_f=u_f[nomatch_from_f_in_i].size;
    u=np.zeros(Nmatched+Nunmatched_i+Nunmatched_f);
    
    # first compute for objects in (i) with a match in (f)
    if (Nmatched>0):
        #u[0:Nmatched] = u_i[sub_id_i] + (u_f[sub_id_f]-u_i[sub_id_i])*tau
        u[0:Nmatched] = u_i[sub_id_i] + (u_f[sub_id_f]-u_i[sub_id_i])*tau
        
        u0 = np.maximum(u_i[sub_id_i],1.e-10)
        u1 = np.maximum(u_f[sub_id_f],1.e-10)
        u_tmp = u0 * np.exp(tau * np.log(u1/u0))
        u[0:Nmatched] = u_tmp 

        if(hsml_rescale==1):
            u0 = 1.*u_i[sub_id_i]; u1 = 1.*u_f[sub_id_f];
            u_tmp = u0*(1.-tau) + u1*tau
            if(pos_interp_flag.size > 1):
                u0=u0[pos_interp_flag]; u1=u1[pos_interp_flag];
                u_tmp[pos_interp_flag] = np.sqrt(u0*u0*(1.-tau) + u1*u1*tau);
            u[0:Nmatched] = u_tmp
        if(mass_rescale==1):
            u_tmp = u_i[sub_id_i] + (u_f[sub_id_f]-u_i[sub_id_i])*tau
            if(pos_interp_flag.size > 1):
                depth = 0.1 # fraction to diminish, which sets width of dip
                #u_tmp_1 = 1.*u_tmp; u_tmp *= 0.;
                u_tmp_1 = 0.*u_tmp; 
                #u_tmp[pos_interp_flag] = u_tmp_1[pos_interp_flag]
                #u_tmp[pos_interp_flag] /= (1. + tau*(1.-tau) * 4.*(1.-depth)/depth)
            u[0:Nmatched] = u_tmp
            
    ## now unmatched: first those in (i) with no match in (f)
    if(Nunmatched_i>0):
        if(mass_rescale==1):
            u_t = 1.*u_i[nomatch_from_i_in_f] * ((1.-tau)**(2.))
        elif(hsml_rescale==1):
            u_t = u_i[nomatch_from_i_in_f] / ((1.-tau)**2.+1.0e-5)
        else:
            u_t = u_i[nomatch_from_i_in_f]
        u[Nmatched:Nmatched+Nunmatched_i] = u_t;
	## now unmatched: now those in (f) with no match in (i)
    if(Nunmatched_f>0):
        if(mass_rescale==1):
            #u_t = u_f[nomatch_from_f_in_i] * (tau**(2.))
            u_t = 1.*u_f[nomatch_from_f_in_i] * (tau**(3.))
        elif(hsml_rescale==1):
            u_t = u_f[nomatch_from_f_in_i] / ((1.-(1.-tau))**2.+1.0e-5)
        else:
            u_t = u_f[nomatch_from_f_in_i]
        u[Nmatched+Nunmatched_i:Nmatched+Nunmatched_i+Nunmatched_f] = u_t;

    gc.collect()
    return u;




########################################################







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
        snapshot_name='snapshot',four_char=four_char)
    if((fname=='NULL')|(fname_ext!='.hdf5')): return dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum;
    gc.collect()

    ## load header information 
    file = h5py.File(fname,'r') # Open hdf5 snapshot file
    header_master = file["Header"] # Load header dictionary (to parse below)
    header_toparse = header_master.attrs
    numfiles = header_toparse["NumFilesPerSnapshot"]
    npartTotal = header_toparse["NumPart_Total"]
    if(np.sum(npartTotal[ptypes])<1): return dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum;
    time = header_toparse["Time"]
    hubble = header_toparse["HubbleParam"]
    if(cosmological==0): 
        ascale = 1.0;
    else:
        ascale = 1.0 * time
        time = cosmological_time(ascale,h=hubble)
    cen_code = np.array(CENTER_BOX) * hubble;
    d_thold = 1.e10
    if(BOX_MAX != 0): d_thold = BOX_MAX * hubble;
    file.close()
    
    ## start parsing the actual snapshot
    xyz=np.zeros((0,3)); xyz_complete=xyz; vxyz=xyz; q=np.zeros((0)); c_all=q; h_all=q; m_all=q; zm_all=q; id_all=np.zeros((0),dtype=long)
    for i_file in range(numfiles):
        if (numfiles>1): fname = fname_base+'.'+str(i_file)+fname_ext  
        if(os.stat(fname).st_size>0):
            file = h5py.File(fname,'r') # Open hdf5 snapshot file
            npart = file["Header"].attrs["NumPart_ThisFile"]
            for ptype in ptypes:
                if(npart[ptype] > 0):
                    p_name_0 = 'PartType'+str(ptype)+'/'
                    xyz_all = np.array(file[p_name_0+'Coordinates'])
                    d = np.amax(np.abs(xyz_all-cen_code),axis=1)
                    ok = np.where(d < d_thold)[0]
                    n0 = ok.shape[0]
                    if(GAS==0):
                        xyz_complete = np.concatenate([xyz_complete,xyz_all])
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
                            dt = np.maximum(dt , min_stellar_age)
                            c_all = np.concatenate([c_all,dt])
                    gc.collect()
                gc.collect()
            file.close()
        gc.collect()
    gc.collect()
    xyz *= ascale / hubble
    m_all /= hubble
    vxyz *= np.sqrt(ascale)
    if(GAS==1):
        h_all *= ascale / hubble
        if(SWAP_TEMP_RHO==1): c_all *= hubble*hubble / (ascale*ascale*ascale)
    else:
        xyz_complete = (xyz_complete-cen_code) * ascale / hubble
        h_all_full = load_allstars_hsml_formovie(snapdir,snapnum,cosmo=cosmological, \
            use_rundir=use_rundir,four_char=four_char,use_h0=h0,
            filename_prefix=filename_prefix,xyzset=xyz_complete);
        d = np.amax(np.abs(xyz_complete),axis=1)
        ok = np.where(d < d_thold  * ascale / hubble)[0]
        h_all = np.concatenate([h_all,h_all_full[ok]])
	## correct to same ID as original gas particle for new stars, if bit-flip applied
    if ((np.min(id_all)<0) | (np.max(id_all)>1.e9)):
        bad = (id_all < 0) | (id_all > 1.e9)
        id_all[bad] += (long(1) << 31)

    ## make sure everything is cast correctly
    ok = (h_all > 0) & (m_all > 0) & (h_all < 100.)
    id_all=id_all[ok]; m_all=m_all[ok]; xyz=xyz[ok,:]; vxyz=vxyz[ok,:]; 
    h_all=1.25*h_all[ok]; c_all=c_all[ok]; zm_all=zm_all[ok];
    gc.collect()
    return id_all, m_all, xyz[:,0], xyz[:,1], xyz[:,2], vxyz[:,0], vxyz[:,1], vxyz[:,2], h_all, c_all, zm_all;


    
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





#########################################################









#
# given a set of positions, project them to the coordinates they 'would be at'
#   with respect to a given camera (for ray-tracing purposes)
#
def coordinates_project_to_camera(x, y, z, \
        camera_pos=[0.,0.,0.], camera_dir=[0.,0.,1.], \
        screen_distance=1.0 ):
    camera_pos=np.array(camera_pos,dtype='d');
    camera_dir=np.array(camera_dir,dtype='d');
    ## recenter at camera location
    x-=camera_pos[0]; y-=camera_pos[1]; z-=camera_pos[2];
    ## line-of-sight distance along camera viewing angle
    znew = x*camera_dir[0] + y*camera_dir[1] + z*camera_dir[2]; 
    ## get perpendicular axes to determine position in 'x-y' space:
    cperp_x,cperp_y = return_perp_vectors(camera_dir);
    xnew = x*cperp_x[0] + y*cperp_x[1] + z*cperp_x[2];
    ynew = x*cperp_y[0] + y*cperp_y[1] + z*cperp_y[2];

    ## project to 'screen' at fixed distance along z-axis :: 
    znew[znew==0.]=-1.0e-10; ## just to prevent nans
    ## need correction factors for other quantities (like hsml, for example)
    r_corr = screen_distance/znew;
    xnew *= r_corr; ## put positions 'at the screen'
    ynew *= r_corr; ## put positions 'at the screen'
    znew = -znew; ## so values in the direction of camera are *negative* (for later sorting)

    gc.collect()
    return xnew, ynew, znew, r_corr


def coordinates_rotate(x_all, y_all, z_all, theta, phi, coordinates_cylindrical=0):
    ## set viewing angle for image
    ## project into plane defined by vectors perpendicular to the vector of this angle
    x = np.cos(phi)*x_all + np.sin(phi)*y_all + 0.*z_all
    y = -np.cos(theta)*np.sin(phi)*x_all + np.cos(theta)*np.cos(phi)*y_all + np.sin(theta)*z_all
    z =  np.sin(theta)*np.sin(phi)*x_all - np.sin(theta)*np.cos(phi)*y_all + np.cos(theta)*z_all
    if (coordinates_cylindrical==1):
        x=x/abs(x)*np.sqrt(x*x+z*z); y=y; z=z
    gc.collect()
    return x,y,z
  

## simple function to wrap x/y/z coordinates near the box edges, to give 
##   the 'correct' x-x0 in periodic coordinates
def periodic_dx( x, x0, box ):
    b0 = box / 2.
    dx = x - x0
    too_large = (dx > b0)
    dx[too_large] -= box
    too_small = (dx < -b0)
    dx[too_small] += box
    gc.collect()
    return dx
    
def get_closest_periodic_value( x_expected, x_final, x_initial, limiter=1.0 ):
    x_remainder = np.mod(x_expected, 2.*np.pi)
    dx = periodic_dx( x_final, x_remainder, 2.*np.pi ) ## returns wrapped dx from x_remainder
    dx_full = x_expected + dx - x_initial
    
    return dx_full + x_initial
    
    # now put a limiter so we don't allow some crazy number of 'laps':
    sign_dx = 0.*dx_full + 1.
    sign_dx[dx_full<0.]=-1.
    dx_full_abs = (sign_dx*dx_full) 
    dx_full_cycles = dx_full_abs / (2.*np.pi) # number of 'laps'
    dx_full_remainder = np.mod(dx_full_cycles, 1.0)
    
    dx_corrected = np.copy(dx_full)
    lminusone = limiter - 1.
    hi = (dx_full_cycles > limiter) & (dx_full_remainder > lminusone)
    dx_corrected[hi] = 2.*np.pi*sign_dx[hi]*dx_full_remainder[hi]
    hi = (dx_full_cycles > limiter) & (dx_full_remainder < lminusone)
    dx_corrected[hi] = 2.*np.pi*sign_dx[hi]*(1.+dx_full_remainder[hi])
    
    gc.collect()
    return dx_corrected + x_initial


def xyz_to_cylindrical(x,y,z,vx,vy,vz):
    for q in [x,y,z,vx,vy,vz]: q=np.array(q,dtype='d')+1.0e-10

    ## get the angular momentum vector (to use for interpolation)
    r = np.sqrt( x*x + y*y + z*z )
    pos = np.zeros([x.size,3],dtype='d'); 
    pos[:,0]=x/r; pos[:,1]=y/r; pos[:,2]=z/r;
    v = np.sqrt( vx*vx + vy*vy + vz*vz )
    vel = np.zeros([vx.size,3],dtype='d'); 
    vel[:,0]=vx/v; vel[:,1]=vy/v; vel[:,2]=vz/v;
    
    j_hat = cross( pos, vel )
    jj = np.sqrt( j_hat[:,0]*j_hat[:,0] + j_hat[:,1]*j_hat[:,1] + j_hat[:,2]*j_hat[:,2] )
    for j in [0,1,2]: j_hat[:,j] /= jj
    
    s = np.argsort(r)
    r_s = r[s]; 
    j_s = np.zeros((r.size,3))
    j_x = np.cumsum(j_hat[s,0] * jj[s])
    j_y = np.cumsum(j_hat[s,1] * jj[s])
    j_z = np.cumsum(j_hat[s,2] * jj[s])
    j_0 = np.sqrt(j_x*j_x + j_y*j_y + j_z*j_z)
    j_x /= j_0; j_y /= j_0; j_z /= j_0;
    j_hat[s,0] = j_x
    j_hat[s,1] = j_y
    j_hat[s,2] = j_z

    ## define the absolute 'phi' relative to an arbitrary (but fixed) 90-degree rotation of the angular momentum vector j
    x_jhat = 0.*j_hat
    x_jhat[:,0]=0.*j_hat[:,0]; x_jhat[:,1]=-j_hat[:,2]; x_jhat[:,2]=j_hat[:,1]
    jj = np.sqrt( x_jhat[:,0]*x_jhat[:,0] + x_jhat[:,1]*x_jhat[:,1] + x_jhat[:,2]*x_jhat[:,2] )
    for j in [0,1,2]: x_jhat[:,j] /= jj
    ## generate y-vector by cross-product a x b = c (gaurantees right-hand rule)
    y_jhat = cross( j_hat, x_jhat )
    jj = np.sqrt( y_jhat[:,0]*y_jhat[:,0] + y_jhat[:,1]*y_jhat[:,1] + y_jhat[:,2]*y_jhat[:,2] )
    for j in [0,1,2]: y_jhat[:,j] /= jj
    
    z_vec = j_hat
    x_vec = x_jhat
    y_vec = y_jhat
    
    xn = x*x_vec[:,0] + y*x_vec[:,1] + z*x_vec[:,2]
    yn = x*y_vec[:,0] + y*y_vec[:,1] + z*y_vec[:,2]
    zn = x*z_vec[:,0] + y*z_vec[:,1] + z*z_vec[:,2]
    vxn = vx*x_vec[:,0] + vy*x_vec[:,1] + vz*x_vec[:,2]
    vyn = vx*y_vec[:,0] + vy*y_vec[:,1] + vz*y_vec[:,2]
    vzn = vx*z_vec[:,0] + vy*z_vec[:,1] + vz*z_vec[:,2]
    
    R = np.sqrt(xn*xn + yn*yn)
    V_R = (vxn * xn + vyn * yn) / R
    V_Phi = (-vxn * yn + vyn * xn) / R
    Phi = np.arctan2( yn , xn )
    ## reset to a 0-2pi system instead of numpy's -pi,pi system
    lo = (Phi < 0.); Phi[lo] += 2.*np.pi

    gc.collect()
    return R, V_R, zn, vzn, Phi, V_Phi, j_hat


    
def xyz_to_polar(x,y,z,vx,vy,vz):
    for q in [x,y,z,vx,vy,vz]: q=np.array(q,dtype='d')+1.0e-10
    x += 1.0e-5; vy += 1.0e-5;

    ## get the angular momentum vector (to use for interpolation)
    r = np.sqrt( x*x + y*y + z*z )
    pos = np.zeros([x.size,3],dtype='d'); 
    pos[:,0]=x/r; pos[:,1]=y/r; pos[:,2]=z/r;
    v = np.sqrt( vx*vx + vy*vy + vz*vz )
    vel = np.zeros([vx.size,3],dtype='d'); 
    vel[:,0]=vx/v; vel[:,1]=vy/v; vel[:,2]=vz/v;
    
    j_hat = cross( pos, vel )
    jj = np.sqrt( j_hat[:,0]*j_hat[:,0] + j_hat[:,1]*j_hat[:,1] + j_hat[:,2]*j_hat[:,2] )
    for j in [0,1,2]: j_hat[:,j] /= jj
    
    ## get spherical polar coordinates
    v_r = vx*pos[:,0] + vy*pos[:,1] + vz*pos[:,2]

    ## now get the vector perpendicular to both j and r
    phi_hat = cross( j_hat, pos )
    jj = np.sqrt( phi_hat[:,0]*phi_hat[:,0] + phi_hat[:,1]*phi_hat[:,1] + phi_hat[:,2]*phi_hat[:,2] )
    for j in [0,1,2]: phi_hat[:,j] /= jj

    v_phi = vx*phi_hat[:,0] + vy*phi_hat[:,1] + vz*phi_hat[:,2]
    v_phi /= r

    ## define the absolute 'phi' relative to an arbitrary (but fixed) 90-degree rotation 
    ##   of the angular momentum vector j
    x_jhat = 0.*j_hat
    x_jhat[:,0]=0.*j_hat[:,0]; x_jhat[:,1]=-j_hat[:,2]; x_jhat[:,2]=j_hat[:,1]
    jj = np.sqrt( x_jhat[:,0]*x_jhat[:,0] + x_jhat[:,1]*x_jhat[:,1] + x_jhat[:,2]*x_jhat[:,2] )
    for j in [0,1,2]: x_jhat[:,j] /= jj
    
    ## generate y-vector by cross-product a x b = c (gaurantees right-hand rule)
    y_jhat = cross( j_hat, x_jhat )
    jj = np.sqrt( y_jhat[:,0]*y_jhat[:,0] + y_jhat[:,1]*y_jhat[:,1] + y_jhat[:,2]*y_jhat[:,2] )
    for j in [0,1,2]: y_jhat[:,j] /= jj

    ## now project r onto this, to obtain the components in this plane
    x_for_phi = x_jhat[:,0]*pos[:,0] + x_jhat[:,1]*pos[:,1] + x_jhat[:,2]*pos[:,2]
    y_for_phi = y_jhat[:,0]*pos[:,0] + y_jhat[:,1]*pos[:,1] + y_jhat[:,2]*pos[:,2]
    phi = np.arctan2( y_for_phi , x_for_phi )
    ## reset to a 0-2pi system instead of numpy's -pi,pi system
    lo = (phi < 0.); phi[lo] += 2.*np.pi

    gc.collect()
    return r, v_r, phi, v_phi, j_hat






###################################################






def test_centering(four_char=0,cosmological=1,
    snapdir='/scratch/01799/phopkins/m12i_hybrid_test/m12_awetzel/m12m/m12m_ref12', #"""location of snapshots"""
    outputdir_master='/work/01799/phopkins/images/', #"""parent folder for frames dump"""
    subfolder_ext='output',force_rebuild=0,snum_min=0,snum_max=1000): #"""subfolder of snapdir containing the snapshots"""

    ## parse the directory names (for file-naming conventions)
    s0=snapdir.split("/"); 
    snapdir_specific=s0[len(s0)-1]; n_s0=1; 
    if(len(snapdir_specific)<=1): snapdir_specific=s0[len(s0)-2]; n_s0=2;
    snapdir_master=''; 
    for j in s0[0:len(s0)+1-n_s0]: snapdir_master += str(j)+'/';
    outputdir_master+='/' ## just to be safe
    full_snapdir = snapdir+'/'+subfolder_ext

    ## to avoid 'camera jitter', the safest thing (albeit more expensive) is to 
    ##   pre-compute the movie center at all snapshots, then use a smoothed tracking of it
    ##   -- call subroutine to check if its done for us, and if not, build it --
    print('... building/loading camera positions list ...')
    i=0; imax=750; snums=[-1];
    while (i < imax):
        fname,fname_base,fname_ext = check_if_filename_exists(full_snapdir,i,four_char=four_char)
        if(fname!='NULL'):
            snums=np.concatenate((snums,np.array([i])));
            if (i > imax-50): imax += 50;
        i += 1;
    snums=snums[snums>=0];
    build_camera_centering(snapshot_list,full_snapdir,snapdir_specific,outputdir_master, 
        force_rebuild=force_rebuild,four_char=four_char,cosmological=cosmological,snum_min=snum_min,snum_max=snum_max)



def get_precalc_zoom_center( snapdir, snapdir_for_naming , output_dir, time_desired, cosmological=0 ):
    fname = get_camera_centering_filename( snapdir_for_naming , output_dir )
    cenfile = open(fname,'r'); brick=np.loadtxt(cenfile).reshape(-1,4)
    time=brick[:,0]; x0=brick[:,1]; y0=brick[:,2]; z0=brick[:,3]; ## define columns
    if (cosmological==1): 
        for q in [x0,y0,z0]: q /= time ## correct to comoving for spline

    nsmoo = 13
    windowfun = 'flat' ## boxcar smoothing
    #nsmoo = 21
    #windowfun = 'hanning' ## gaussian-like smoothing (should use larger window function)
    cen = np.zeros(3)
    for j,p in zip([0,1,2],[x0,y0,z0]):
        pp = util.smooth(p, window_len=nsmoo, window=windowfun)
        #pp = p ## linear interpolation (in case we don't have enough points) #DEBUG!
        if (cosmological==1): 
            pp_interp = time_desired * np.interp( np.log(time_desired), np.log(time), pp )
        else:
            pp_interp = np.interp( time_desired, time, pp )
        cen[j] = pp_interp

    print('... ... ... getting center for frame at t=',time_desired,' : ',cen)
    return cen
 


def build_camera_centering( snums, snapdir , snapdir_for_naming, output_dir, 
    center_on_bh = 0, # center on first BH particle (good for nuclear sims)
    center_on_com = 0, # center on center of mass (good for merger sims)
    set_fixed_center = 0, # option to set a fixed center
    force_rebuild = 0, # force the routine to build from scratch (don't check for file)
    use_h0=1,four_char=0,cosmological=1,skip_bh=0,
    min_searchsize=1.,search_deviation_tolerance=0.03,ptype=4,
    snum_min=0,snum_max=1000):

    ## first check if the tabulated centers file already exists:
    fname = get_camera_centering_filename(snapdir_for_naming , output_dir )
    if(force_rebuild==0): 
        if(os.path.exists(fname)): return 1

    snums=snums[(snums>=snum_min)&(snums<=snum_max)]
    ## ok, if we're here, center-table doesn't exist, so let's build it:
    n_snums = np.array(snums).size
    time_all = np.zeros( (n_snums) )
    cen = np.zeros( (n_snums,3) )
    outfi_cen = open(fname,'a')
    set_fixed_center = np.array(set_fixed_center)
    center_on_zero = 0
    n_prev = 0
    if (set_fixed_center.size==1): 
        if (set_fixed_center==1): center_on_zero=1
    for i in range(n_snums):
        print('centering for snap ',snums[i],' in ',snapdir)
        
        fname,fname_base,fname_ext = check_if_filename_exists(snapdir,snums[i],four_char=four_char)
        file = h5py.File(fname,'r')
        time_all[i] = file["Header"].attrs["Time"]
        hubble = file["Header"].attrs["HubbleParam"]
        file.close()
        cen[i,:] = [0.,0.,0.]

        if (center_on_zero==1): 
            continue
        if (set_fixed_center.size==3):
            cen[i,:]=set_fixed_center
            continue

        max_searchsize = 1000.
        cen_guess=np.zeros(3)
        rcorr = 0.*time_all + hubble
        if(cosmological==1): rcorr /= time_all
        
        if(i > 0):
            cen_guess = cen[i-1,:] * rcorr[i-1]
            if((i > 1)&(n_prev>100.)):
                d_cen = np.sqrt(np.sum((cen[i-1,:]*rcorr[i-1]-cen[i-2,:]*rcorr[i-2])**2))
                #thold = 8.*d_cen/rcorr[i]
                thold = 5.*d_cen/rcorr[i]
                if(thold < max_searchsize): max_searchsize=thold
                #if(max_searchsize < 15.*min_searchsize  ): max_searchsize=15.*min_searchsize
                if(max_searchsize < 10.*min_searchsize  ): max_searchsize=10.*min_searchsize

        ## if not using any of the methods above, use our fancier iterative 
        ##   centering solution for the 'dominant' structure
        cen[i,:], n_prev = fast_zoom_center( snapdir, snums[i], \
            max_searchsize=max_searchsize, min_searchsize=min_searchsize, search_deviation_tolerance=search_deviation_tolerance,
            cen_guess=cen_guess,ptype=ptype,four_char=four_char,cosmological=cosmological)
        print('..zoom_cen at ',cen[i,:])
        sys.stdout.flush()
        gc.collect()
        v = np.array([ time_all[i], cen[i,0], cen[i,1], cen[i,2] ]).reshape(-1,4)
        np.savetxt(outfi_cen,v)
    gc.collect()
    ## all done, write out these results
    outfi_cen.close()
    return 1



def fast_zoom_center(sdir, snum, snapshot_name='snapshot', extension='.hdf5', four_char=0,
    max_searchsize=1000., min_searchsize=1., search_deviation_tolerance=0.05,
    cen_guess=[0.,0.,0.], ptype=4, cosmological=1):
    
    fname,fname_base,fname_ext = check_if_filename_exists(sdir,snum,\
        snapshot_name=snapshot_name,extension=extension,four_char=four_char)
    cen = np.array(cen_guess)
    if(fname=='NULL'): return cen
    if(fname_ext!='.hdf5'): return cen
    gc.collect()
    file = h5py.File(fname,'r') # Open hdf5 snapshot file
    header_master = file["Header"] # Load header dictionary (to parse below)
    header_toparse = header_master.attrs
    numfiles = header_toparse["NumFilesPerSnapshot"]
    npartTotal = header_toparse["NumPart_Total"]
    if(npartTotal[ptype]<1): return cen, 0    
    if(npartTotal[ptype]<1): ptype = 0
    if(npartTotal[ptype]<1): ptype = 1
    boxsize = header_toparse["BoxSize"]
    time = header_toparse["Time"]
    if(cosmological==0): time=1
    hubble = header_toparse["HubbleParam"]
    file.close()
    rcorr = hubble / time; # scale to code units
    max_searchsize *= rcorr; min_searchsize *= rcorr; search_deviation_tolerance *= rcorr;
    if(max_searchsize > boxsize): max_searchsize=boxsize
    d_thold = max_searchsize
    if(cen[0]==0.): d_thold=1.e10
    gc.collect()
    niter = 0
    xyz = np.zeros((0,3))
    while(1):
        xyz_m=np.zeros(3); n_t=np.zeros(1);
        if(niter == 0):
            for i_file in range(numfiles):
                if (numfiles>1): fname = fname_base+'.'+str(i_file)+fname_ext  
                if(os.stat(fname).st_size>0):
                    file = h5py.File(fname,'r') # Open hdf5 snapshot file
                    npart = file["Header"].attrs["NumPart_ThisFile"]
                    if(npart[ptype] > 1):
                        xyz_all = np.array(file['PartType'+str(ptype)+'/Coordinates/'])
                        d = np.amax(np.abs(xyz_all-cen),axis=1)
                        ok = np.where(d < d_thold)
                        n0 = ok[0].shape[0]
                        xyz_all = xyz_all.take(ok[0],axis=0)
                        if(xyz_all.size > 0): 
                            xyz = np.concatenate([xyz,xyz_all])
                            xyz_m += np.sum(xyz_all,axis=0)
                            n_t += n0
                        gc.collect()
                    file.close()
                gc.collect()
            gc.collect()
            xyz_prev = xyz
        else:
            d = np.amax(np.abs(xyz_prev-cen),axis=1)
            ok = np.where(d < d_thold)
            n0 = ok[0].shape[0]
            xyz = xyz_prev.take(ok[0],axis=0)
            if(n0 > 1):
                xyz_m += np.sum(xyz,axis=0)
                n_t += n0
        gc.collect()
        niter += 1;
        if(n_t[0] <= 0): 
            d_thold *= 1.5
        else:
            xyz_m /= n_t[0]
            d_cen = np.sqrt(np.sum((cen-xyz_m)**2))
            print('cen_o=',cen[0],cen[1],cen[2],' cen_n=',xyz_m[0],xyz_m[1],xyz_m[2],' in box=',d_thold,' cen_diff=',d_cen,' min_search/dev_tol=',min_searchsize,search_deviation_tolerance)
            if(niter > 100): break
            if(d_thold <= min_searchsize): break
            if(d_cen <= search_deviation_tolerance): break
            if(n_t[0] <= 10): break
            cen = xyz_m
            d_thold /= 5.1
            if(d_thold < 2.*d_cen): d_thold = 2.*d_cen
            if(max_searchsize < d_thold): d_thold=max_searchsize
            xyz_prev = xyz
        gc.collect()
    gc.collect()
    return xyz_m / rcorr, npartTotal[ptype]



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
        if(os.stat(fname).st_size <= 0):
            ## file exists but is null size, do not use
            fname_found = 'NULL'
            fname_base_found = 'NULL'
            fname_ext = 'NULL'
            continue;
        fname_found = fname;
        fname_base_found = fname_base;
        fname_ext = extension_touse
        break; # filename does exist! 
    return fname_found, fname_base_found, fname_ext;



def get_camera_centering_filename(snapdir_for_naming , output_dir):
    ext = '.mv_cen'
    fname = output_dir + snapdir_for_naming + ext
    return fname



def get_snapshot_precomputed_z0_cen(sdir):
    cen=[0.,0.,0.]
    if(sdir=='m09_ref11_adapt'): cen=[ 4710.58752192,  3425.17425247,  3176.85378332]
    if(sdir=='m09_ref11_adapt_metaldiffmhd_hiThold'): cen=[ 4709.3307902,   3425.09992071,  3178.07732329]
    if(sdir=='m09_ref11_fixedsoft'): cen=[ 4707.89059383, 3424.72304372, 3178.41580743]
    if(sdir=='sf_model_tests_default'): cen=[ 41828.52813216,  44147.23738634,  46238.4724623 ]
    if(sdir=='m12i_ref11'): cen=[ 41718.11654436,  44229.69343256,  46423.92902392]
    if(sdir=='m12i_ref11_tsph'): cen=[ 41817.80042241,  44157.58171149,  46286.9059771 ]
    if(sdir=='m10v_ref8'):  cen=[ 3403.30284402,  3298.51326896,  3529.32402938]
    if(sdir=='m10v_ref9'):  cen=[ 3397.61109927,  3290.66750332,  3524.44558671]
    if(sdir=='m10v_ref10'): cen=[ 3418.28945977,  3291.88373753,  3566.75524528]
    if(sdir=='m10v_ref11_adapt'): cen=[ 3427.06370651,  3271.65994972,  3523.04310327]
    if(sdir=='m10v_ref11_adapt_metaldiffmhd'): cen=[ 3426.16529536,  3271.85564295,  3522.90864404]
    if(sdir=='m10v_ref11'): cen=[ 3427.06370651,  3271.65994972,  3523.04310327]
    if(sdir=='m10q_ref8'):  cen=[ 3335.72293591,  2974.47073640,  3350.93201216]
    if(sdir=='m10q_ref9'):  cen=[ 3329.01555615,  2968.60596968,  3343.08745896]
    if(sdir=='m10q_ref10'): cen=[ 3329.05788868,  2985.49113833,  3311.77384541]
    if(sdir=='m10q_ref11'): cen=[ 3313.02295174,  2945.93301173,  3394.93218315]
    if(sdir=='m10q_ref11_fixed'): cen=[ 3316.39482422, 2943.0249495,  3393.59247881]
    if(sdir=='m10q_ref11_adapt_alldiffmhd'): cen=[ 3303.94783806,  2947.05529221,  3397.22271046]
    if(sdir=='m10q_ref11_adapt_hiThold'): cen=[ 3313.13005136,  2943.92293291,  3397.17398351]
    if(sdir=='m12i_ref11'): cen=[ 41718.11654436,  44229.69343256,  46423.92902392]
    if(sdir=='m12i_ref12'): cen=[ 41828.29846686,  44131.82912578,  46298.30386540]
    if(sdir=='m12i_ref13'): cen=[ 41875.75676213,  44122.38220725,  46257.48356735]
    if(sdir=='m12i_ref12_adapt_metaldiff2'): cen=[41830.50437144,  44239.39314392,  46234.6609263]
    if(sdir=='m12i_ref11_noRad_tsph'): cen=[ 41818.87257281,  44160.69550523,  46289.22666334]
    if(sdir=='m12f_ref12'): cen=[ 38748.68542788,  47643.97050262,  46763.73983554]
    if(sdir=='m12f_ref12_fb-sym'): cen=[ 38748.69056642,  47644.00698554,  46763.80678027]
    if(sdir=='m12f_ref12_fb-iso'): cen=[ 38680.5433167,   47707.30194236,  46818.98373489]
    if(sdir=='m12b_ref12_fb-sym'): cen=[ 39311.00446494, 41979.01997792, 39152.44726635]
    if(sdir=='m12b_ref12_fb-iso'): cen=[ 39337.51460655,  41756.55378332,  39179.1354953 ]
    if(sdir=='m12c_ref12_fb-sym'): cen=[ 36008.06717868,  49152.75510195,  46821.18188375]
    if(sdir=='m12c_ref12_fb-iso'): cen=[ 35916.69167078,  49290.06782276,  46877.39851561]
    if(sdir=='m12m_ref12_fb-sym'): cen=[ 37561.37411999,  41583.64931681,  46769.30807611]
    if(sdir=='m12m_ref12_fb-iso'): cen=[ 37553.97725373,  41653.68221635,  46756.85800067]
    if(sdir=='m12f_ref13_fb-sym'): cen=[ 38744.24350156,  47701.65789227,  46790.76338099]
    if(sdir=='m12f_ref13'): cen=[ 38744.24350156,  47701.65789227,  46790.76338099]
    if(sdir=='m11v_ref13_adapt'): cen=[ 42824.85443495,  44154.3228809,   43528.57110144]
    if(sdir=='m11q_ref13_adapt'): cen=[ 42242.33758432,  45740.82325377,  40512.26357658]
    if(sdir=='m12i_ref12_fire1'): cen=[ 41925.08412405,  44677.04569901,  46068.99276398]
    if(sdir=='m09_ref11_fire1'): cen=[ 4705.95809371,  3424.56268181,  3178.0762833 ]
    if(sdir=='m10q_ref11_fire1'): cen=[ 3207.62480088,  2976.60457247,  3286.54402602]
    if(sdir=='m10v_ref10_fire1'): cen=[ 4494.65998758,  3167.99864101,  2686.66612086]
    if(sdir=='m11v_ref12_fire1'): cen=[ 43229.3474205,   44158.01019037,  43678.13068558]
    if(sdir=='m11q_ref13_fire1'): cen=[ 42051.52966859,  45901.56828362,  39757.28554621]
    if(sdir=='m12i_ref12_fb_iso'): cen=[ 41843.52124705,  44171.67922194,  46253.53395134]
    if(sdir=='m12i_ref11_fb-aniso'): cen=[ 41710.99239422,  44338.27894807,  46346.86432687]
    if(sdir=='m12i_ref11_fb_iso'): cen=[ 41797.51755388,  44322.46511334,  46321.86212181]
    if(sdir=='m12i_ref12_fb_aniso'): cen=[ 41867.35752106,  44164.56941237,  46238.50485002]
    if(sdir=='m10q_ref10_fixedsoft_20'): cen=[ 3330.1619776,   2981.91197077,  3312.14271496]
    if(sdir=='m10q_ref10_fixedsoft_dm'): cen=[ 3318.74676144,  2982.3481015,   3316.68912189]    
    if(sdir=='m10q_ref11_fixedsoft'): cen=[ 3316.39482422,  2943.0249495,   3393.59247881]
    if(sdir=='m10q_ref11_dm_fixedsoft'): cen=[ 3302.74423667,  2947.34243833,  3400.34404091]
    if(sdir=='m10v_ref10_fixedsoft_dm'): cen=[ 3407.94699333,  3299.88783958, 3559.30869286]
    if(sdir=='m10v_ref11_fixedsoft'): cen=[ 3426.82879833,  3271.68589105,  3522.73985386]
    if(sdir=='m10v_ref11_dm_fixedsoft'): cen=[ 3415.57177299,  3277.97178013,  3517.26401981]
    if(sdir=='m12i_ref12_dm'): cen=[ 41824.99971528,  44152.94674268,  46280.59357175]
    if(sdir=='m12i_ref13_dm'): cen=[ 41820.51225595,  44153.30746753,  46272.85728992]
    if(sdir=='m12f_ref11_fb-sym'): cen=[ 38616.98367136,  47776.05817846,  46847.70543562]
    if(sdir=='m12f_ref12_dm'): cen=[ 38723.17211486,  47653.0094322,   46787.35209018]
    if(sdir=='m12f_ref12_dm_rad6'): cen=[ 38724.56059854,  47658.2378309,   46786.7506989 ]
    if(sdir=='m12f_ref13_dm'): cen=[ 38715.93421233,  47652.45341032,  46782.89188415]
    if(sdir=='m11q_ref13_fixedsoft'): cen=[ 42259.63436309,  45737.67365735,  40502.37933148]
    if(sdir=='m11v_ref13_fixedsoft'): cen=[ 42832.21307527,  44182.01152496,  43521.16148749]
    if(sdir=='m12f_ref11_dm_fixedsoft'): cen=[ 38732.23485065,  47642.62676188,  46775.60105084]
    if(sdir=='m12i_ref11_dm_smallsoft'): cen=[ 41827.37934501,  44152.90244138,  46274.98437919]
    if(sdir=='m12i_ref11_dm_small2soft'): cen=[ 41827.73049651,  44151.33876179,  46274.76064483]
    if(sdir=='m12i_ref11_dm'): cen=[ 41827.3526075,   44151.2047719,   46274.29765069]
    if(sdir=='m12i_ref11_dm_bigsoft'): cen=[ 41828.11838077,  44150.94670074,  46274.22638959]
    if(sdir=='m12i_ref11_dm_comoving600pc'): cen=[ 41826.84372491,  44152.31257079,  46274.1704996 ]
    if(sdir=='m12i_ref11_dm_ags'): cen=[ 41828.41312479,  44152.51080231,  46272.72380872]
    if(sdir=='m12i_ref13_dm_res-20pc'): cen=np.array([ 41820.82506815,  44152.98201607, 46272.80317769])
    if(sdir=='m12i_ref13_dm_res-80pc'): cen=np.array([ 41820.71264973,  44153.48068267,  46272.7722742 ])
    if(sdir=='m12i_ref13_dm_res-160pc'): cen=np.array([ 41820.83440886,  44153.3748299,   46272.86645861])
    if(sdir=='m10v_ref8_dm'): cen=[ 3399.08016416,  3303.44614895,  3522.91628842]
    if(sdir=='m10v_ref8_fixedsoft_dm'): cen=[ 3396.16711927,  3301.70932492,  3525.64239658]
    if(sdir=='m10v_ref9_dm'): cen=[ 3388.88283372,  3297.18020894,  3517.89061786]
    if(sdir=='m10v_ref9_fixedsoft_dm'): cen=[ 3388.85220596,  3297.1886898,   3517.74826294]
    if(sdir=='m10v_ref10_dm'): cen=[ 3408.40547067,  3299.4954469,   3558.31712581]
    if(sdir=='m10v_ref10_fixedsoft_dm'): cen=[ 3407.94699333,  3299.88783958,  3559.30869286]
    if(sdir=='m10v_ref11_dm_adapt'): cen=[ 3415.89241673,  3277.89113139,  3516.8104357 ]
    if(sdir=='m12i_ref13_fb_aniso'): cen=[ 41875.12718359,  44112.15174466,  46273.36882759]
    if(sdir=='m12i_ref13_MHD_CV_TD'): cen=[41815.491,44136.821,46267.662]
    return cen







###################################################################


