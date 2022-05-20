import numpy as np
import pfh_utils as util
import os.path
import visualization.image_maker as image_maker
from gizmopy.quicklook import estimate_zoom_center
from gizmopy.load_fire_snap import load_fire_snap, evaluate_if_cosmological
from gizmopy.load_from_snapshot import load_from_snapshot
import sys
# from scipy.interpolate import interp1d

## here is the main routine for the galaxy simulation movies: 
##   it builds the snapshot list, does the centering, goes through frames and calls
##   the appropriate routines for interpolation between snapshots, and then calls the 
##   image_maker routine to render the frame: many switches here (most don't need to be set), 
##   but tune it as you need to (I recommend a low snapshot spacing to test)
##
def movie_maker(\
    sdir, #='output', #"""location of snapshots"""
    output_directory, #='movie', #"""parent folder for frames dump"""
    field_of_view=30., #"""scale (in code units) of movie box size"""
    center_filename='movie_center.txt', # name of the file that contains scalefactor vs time in snapshot directory
    image_key='star', #"""determines if image is 'stars','gas','xr' etc, and overlay layers"""
    frames_per_gyr=100., #"""sets spacing of frames (even in time)"""
    use_physical_units=True, #"""convert to physical units?"""
    time_min=0, #"""initial time of movie"""
    time_max=0, #"""final time (if<=0, just set to time of maximum snapshot)"""
    frame_min=0, #"""initial frame number to process"""
    frame_max=1.0e10, #"""final frame number to process"""
    snap_min=0, #"""minimum snapshot number in list to scan"""
    snap_max=0, #"""maximum snapshot number to scan (if =0, will just use largest in snapdir)"""
    convert_float32=True, #convert all properties to float32
    saturation_level=None,  # saturation_level to pass to image_maker
    overwrite=False,
    **kwargs,  # passed to image_maker, parse_gasvals_bands_toplot, load_from_snapshot, etc.
    ):
    """
    make a movie of a simulation.

    calculates the center at each time, loads snapshots in i, i+1 pairs and 
    interpolates between, using a smoothed center from the pre-calculated track, 
    to make a movie that's smooth in time.  calls image_maker.image_maker on 
    the interpolated particle positions.

    Args:
        sdir (string):  directory containing the snapshot files/snapshot 
                        directories.  typically 'output' or 
                        '<simulation directory>/output'

        output_directory (string):  directory to save the images (and raw 
                                    image dumps) to

        field_of_view (float):  size of the box in physical code units

        center_filename (string):  path to file that contains the scale factor 
                                   and x-y-z position of the center in physical 
                                   code units (will be created using default
                                   centering algorithm if it doesn't exist)

        image_key (string):  type of images to make.  see 
                             image_maker.image_maker for an explanation

        frames_per_gyr (int):  number of frames to make per gigayear

        use_physical_units (bool):  whether or not to make the the 
                                    image in physical (True) or comoving
                                    (False).  makes the biggest difference
                                    in gas movies and large-scale movies, 
                                    where you can see the expansion of the
                                    universe

        time_min, time_max (float):  min and max time (non-cosmological)
                                     or scale factors to make the movie
                                     between.  if max is set to zero, 
                                     then it's ignored.

        frame_min, frame_max (float):  min and max frame numbers to visualize.
                                       note that these are specific to a given
                                       simulation and frames_per_gyr value, so 
                                       be careful with this one.

        snap_min, snap_max (float):  min and max snapshot numbers to include in
                                     the movie.  if max is zero, again ignored.

        convert_float32 (bool):  whether or not to convert the particle data 
                                 to 32-bit to save memory

        saturation_level (float):  saturation level to pass to image_maker.  note
                                   that this should be set positive, since you can't
                                   rely on the dynamically-set saturation levels to be
                                   smooth, so you get weird pulsing in the movies.  
                                   could make this a smooth function in principle,
                                   but that's not supported at the moment.  units 
                                   should be code units; for default FIRE, that means 
                                   solar luminosities/kpc^2.  raise this if you 
                                   entire image is saturated, and lower it if your
                                   entire image is dark.

        overwrite (bool):  whether or not to overwrite existing image files, or skip them
        
        kwargs:  any keyword arguments, mostly passed along to image_maker.image_maker
                 (but also sent to load_fire_snap, parse_gasvals_bands_toplot, and a
                 couple others.)  important ones include:
                    * pixels:  sets the size/resolution of the image
                    * dyanmic_range:  sets the dynamic range of the images 
                    * DesNgb:  sets the desired number of neighbors for smoothing the star particles
                    * projection: set to either a blank string or 'camera_tocenter'.  
                                  don't think I recommend using 'camera' without 'tocenter', 
                                  since tht places the camera at the center of the galaxy at each time.
                                  could be interesting though!
                    * camera_opening_angle:  sets the opening angle of the camera lense, if using 
                                             a camera projection 

    Returns:
        nothing, but it saves the images in the output_directory

    """
    # get the default values of arguments to image_maker algorithmically:
    image_maker_defaults = util.get_default_args(image_maker.image_maker)

    ## first decide whether its a gas or stellar image and what gas value we're making
    bands = kwargs.get('bands', image_maker_defaults['bands'])
    kernel_width = kwargs.get('kernel_width', image_maker_defaults['kernel_width'])

    image_key, effects_key, gasvalue_to_image, bands, kernel_width = image_maker.parse_gasvals_bands_toplot(
        image_key, bands, kernel_width, **kwargs)

    stars_are_primary_type = 'star' in image_key
    
    ## set the saturation according to the box size, if not set explicitely:
    if isinstance(saturation_level, str) or saturation_level is None:
        saturation_level = set_saturation_from_fov(field_of_view)
    if saturation_level < 0:
        print("\n-- warning -- saturation_level is set dynamically, which will cause flickering in the final movie")
        print("-- recommend setting by hand, or use None (or any string) to scale to the box size")

    ## build the snapshot list in the directory of interest
    print('... building snapshot list in directory ...')
    snapshot_list = build_snapshot_list(sdir, **kwargs)
    if(snap_max<=0): 
        snap_max=np.max(snapshot_list)

    assert len(snapshot_list), "Didn't find any snapshots in {}".format(sdir)
    cosmological = evaluate_if_cosmological(sdir, snapshot_list[0])

    ## to avoid 'camera jitter', the safest thing (albeit more expensive) is to 
    ##   pre-compute the movie center at all snapshots, then use a smoothed tracking of it
    ##   -- call subroutine to check if its done for us, and if not, build it --
    print('... ensuring camera position exists ...')
    ensure_center_of_time_exists(snapshot_list, sdir, center_filename, **kwargs)

    print('... building list of times for movie frames ...')
    ## now build the times for each frame of the movie
    time_frame_grid, a_scale_grid = build_time_frame_grid(sdir, snapshot_list, 
          frames_per_gyr, time_min=time_min, time_max=time_max, cosmological=cosmological)


    ## now enter main loop over snapshots in the simulation set
    snaps_to_do = (snapshot_list >= snap_min) & (snapshot_list <= snap_max)
    ii = np.arange(snapshot_list.size)
    snap_grid_to_loop = ii[snaps_to_do]
    frame_number_grid = np.arange(time_frame_grid.size)

    print('... entering main loop over {} snapshots ...'.format(snap_grid_to_loop.size))
    for i_snap in snap_grid_to_loop[:-1]:
        ## load header info for snapshot 'i' and 'f = i + 1'
        if cosmological:
            a_scale_i = load_from_snapshot('Time', -1, sdir, snapshot_list[i_snap], **kwargs)
            t_i = cosmological_time(a_scale_i)

            a_scale_f = load_from_snapshot('Time', -1, sdir, snapshot_list[i_snap+1], **kwargs)
            t_f = cosmological_time(a_scale_f)

        else:
            t_i = a_scale_i = load_from_snapshot('Time', -1, sdir, snapshot_list[i_snap], **kwargs)
            t_f = a_scale_i = load_from_snapshot('Time', -1, sdir, snapshot_list[i_snap+1], **kwargs)
            
        ## now calculate whether there are any frames in the time-range between the snapshots
        snapshot_delta_t = t_f - t_i
        print('Timestep from snapshot {} (ti = {:.1f}) - {} ({:.1f}) is {:.2f}'.format(
            snapshot_list[i_snap], t_i, snapshot_list[i_snap+1], t_f, snapshot_delta_t))

        check_frames_to_do = (time_frame_grid >= t_i) & (time_frame_grid <= t_f) & \
            (frame_number_grid <= frame_max) & (frame_number_grid >= frame_min)
        frames_to_do = frame_number_grid[check_frames_to_do]

        ## now only keep the frames where no file yet exists (if we're not overwriting)
        if overwrite==False:
            ## be verbose at this point, since it's a little surprising we got here
            frames_to_do = mask_existing_frames(frames_to_do, sdir, 
                output_directory, field_of_view, image_key, verbose=True, **kwargs)

        ## now we're ready to loop through these frames
        ## note that above check means we skip any snapshots where
        ## all the frames between them are already done!
        n_frames_to_do = frames_to_do.size

        if(n_frames_to_do>0):
            print('... processing {} frames in this snapshot interval ...'.format(
                n_frames_to_do))

            ### get snapshot centers in physical coordinates
            center_i = get_precalc_zoom_center(center_filename, a_scale_i, cosmological=cosmological)
            center_f = get_precalc_zoom_center(center_filename, a_scale_f, cosmological=cosmological)

            print('... ... loading snapshot bricks ... ...')
            ## tolerance for keeping particles outside the box, for this calculation
            xmax_tmp = (field_of_view * 1.5) + 10.

            ## load the actual --data-- for snapshot (i)
            ##  -- will load in physical coordinates, so pas sin the center in physical
            id_i,m_i,x_i,y_i,z_i,vx_i,vy_i,vz_i,h_i,c_i,zm_i = load_pos_vel_etc(
                sdir, snapshot_list[i_snap], CENTER_BOX=center_i, 
                GAS=(stars_are_primary_type==False), BOX_MAX=xmax_tmp, 
                convert_float32=convert_float32, **kwargs)

            if id_i.size <= 3:
                print("not enough primary particles in snapshot {} to image; continuing to next chunk".format(snapshot_list[i_snap]))
                continue

            if stars_are_primary_type:
                id_i_g,m_i_g,x_i_g,y_i_g,z_i_g,vx_i_g,vy_i_g,vz_i_g,h_i_g,c_i_g,zm_i_g = load_pos_vel_etc(
                    sdir, snapshot_list[i_snap], CENTER_BOX=center_i,
                    GAS=True,SWAP_TEMP_RHO=True,BOX_MAX=xmax_tmp, 
                    convert_float32=convert_float32, **kwargs)

            ## correct to comoving coordinates for interpolation (so correctly capture hubble flow)
            ## (splitting initial and final rescalings so no worry of doing it twice on 'recycling' step)
            if cosmological:
                center_i /= a_scale_i

                for vec in [x_i,y_i,z_i,h_i,vx_i,vy_i,vz_i]: 
                    vec /= a_scale_i
                
                if stars_are_primary_type:
                    for vec in [x_i_g,y_i_g,z_i_g,h_i_g,vx_i_g,vy_i_g,vz_i_g]: 
                        vec /= a_scale_i

            ## now load snapshot (f) [different center etc., have to load fresh]
            # again, loads physical coords, so pass center in physical coords too
            id_f,m_f,x_f,y_f,z_f,vx_f,vy_f,vz_f,h_f,c_f,zm_f = load_pos_vel_etc(
                sdir, snapshot_list[i_snap+1], CENTER_BOX=center_f,
                GAS=(stars_are_primary_type==False), BOX_MAX=xmax_tmp, 
                convert_float32=convert_float32, **kwargs)
            
            if id_f.size <= 3:
                print("not enough primary particles in snapshot {} to image; continuing to next chunk".format(snapshot_list[i_snap+1]))
                continue

            if stars_are_primary_type:
                id_f_g,m_f_g,x_f_g,y_f_g,z_f_g,vx_f_g,vy_f_g,vz_f_g,h_f_g,c_f_g,zm_f_g = load_pos_vel_etc(
                    sdir, snapshot_list[i_snap+1], CENTER_BOX=center_f,
                    GAS=True,SWAP_TEMP_RHO=True,BOX_MAX=xmax_tmp,
                    convert_float32=convert_float32, **kwargs)

            ## correct to comoving coordinates for interpolation (so correctly capture hubble flow)
            ## do this regardless of whether or not we're going to plot physical units; that'll be handled later
            if cosmological:
                center_f /= a_scale_f

                for vec in [x_f,y_f,z_f,h_f,vx_f,vy_f,vz_f]: 
                    vec /= a_scale_f
                
                if stars_are_primary_type:
                    for vec in [x_f_g,y_f_g,z_f_g,h_f_g,vx_f_g,vy_f_g,vz_f_g]: 
                        vec /= a_scale_f

            ## before going further, check that we actually have particles to process, 
            ##   and if not skip this frame
            if (m_i.size+m_f.size < 5):
                continue;
        
            print('... ... matching ids in the snapshots ... ...')
            ## predefine matches to save time in the loop below
            sub_id_i,sub_id_f,nomatch_from_i_in_f,nomatch_from_f_in_i = compile_matched_ids(id_i, id_f)
            if stars_are_primary_type:
                sub_id_i_g,sub_id_f_g,nomatch_from_i_in_f_g,nomatch_from_f_in_i_g = compile_matched_ids(id_i_g, id_f_g)

            print('... ... entering frame loop ... ...')
            ## loop over frames in between the two snapshots just pulled up
            for i_frame in range(n_frames_to_do):
                j_of_frame = frames_to_do[i_frame];
                dt = time_frame_grid[j_of_frame] - t_i
                time_of_frame = time_frame_grid[j_of_frame]

                print('... ... ... interpolating for frame ',i_frame+1,'/',n_frames_to_do,'... ... ...')
                ## this is the interpolation step between the two snapshots for each frame
                x_all,y_all,z_all,m_all,h_all,c_all,zm_all = interpolation_for_movie(dt, snapshot_delta_t, 
                    m_i,h_i,c_i,zm_i,m_f,h_f,c_f,zm_f, 
                    x_i,y_i,z_i,vx_i,vy_i,vz_i, x_f,y_f,z_f,vx_f,vy_f,vz_f, 
                    sub_id_i,sub_id_f, nomatch_from_i_in_f,nomatch_from_f_in_i, 
                    center_i,center_f,use_polar_interpolation=0,stars=stars_are_primary_type)

                if stars_are_primary_type:
                    # now interpolate the gas
                    x_all_g,y_all_g,z_all_g,m_all_g,h_all_g,c_all_g,zm_all_g = interpolation_for_movie(dt, snapshot_delta_t, 
                        m_i_g,h_i_g,c_i_g,zm_i_g,m_f_g,h_f_g,c_f_g,zm_f_g, 
                        x_i_g,y_i_g,z_i_g,vx_i_g,vy_i_g,vz_i_g, x_f_g,y_f_g,z_f_g,vx_f_g,vy_f_g,vz_f_g, 
                        sub_id_i_g,sub_id_f_g, nomatch_from_i_in_f_g,nomatch_from_f_in_i_g, 
                        center_i,center_f,use_polar_interpolation=0,stars=0)
                ## ok, now have the interpolated gas(+stellar) quantities needed for each image

                ## transform back to physical coordinates for plotting (if requested)
                a = a_scale_grid[j_of_frame]
                if cosmological and use_physical_units:
                    time_of_frame = a
                    for vec in [x_all,y_all,z_all,h_all]: 
                        vec *= a
                    if stars_are_primary_type:
                        for vec in [x_all_g,y_all_g,z_all_g,h_all_g]: 
                            vec *= a
                ### so now all of our coordinates are in physical coordinates unless use_physical_units==False

                print('... ... ... centering and setting passed variables ... ... ...')
                
                ## interpolate the center between the snapshots in physical coordinates
                cen = get_precalc_zoom_center(center_filename, a_scale_grid[j_of_frame], cosmological=cosmological)
                
                ## if we're not doing physical units, then we transform to comoving
                if cosmological and (use_physical_units==False):  cen = cen / a
                ## so now my center is also in physical units unless use_physical_units==False

                ## now stuff interpolated data into a dictionary for the image maker to pop out
                ### remember -- don't pre-process the data!  let image_maker handle all that
                #### slight correction -- we'll recenter, but not do anything else
                particle_data = {}

                # store the time in there as a scalar
                particle_data['Time'] = a_scale_grid[j_of_frame]

                # first, positions -- main particle type is always just xyz
                ## go ahead and shift to have the center at 0, 0, 0, and pass cen as 0, 0, 0
                particle_data['xyz'] = np.vstack((x_all,y_all,z_all)).T - cen

                # if we have stars, then the gas is xyz_g
                if stars_are_primary_type:
                    particle_data['xyz_g'] = np.vstack((x_all_g,y_all_g,z_all_g)).T - cen

                # now do everything else I need, starting with the main particle type
                # if that's gas, then I prefix it with gas_
                prefix = 'gas_'
                if stars_are_primary_type:
                    # otherwise, no prefix
                    prefix = ''

                particle_data[prefix+'SmoothingLength'] = h_all
                particle_data[prefix+'Masses'] = m_all
                particle_data[prefix+'Z'] = zm_all

                # now handle the "color" datasets:
                if stars_are_primary_type:
                    # if stars, then c_all is stellar age and c_all_g is gas temperature
                    particle_data['StellarAgeGyr'] = c_all
                else:
                    # if gas, then c_all is gas temperature (or whatever) and c_all_g is undefined
                    particle_data['gas_'+gasvalue_to_image] = c_all

                # now handle the gas if I'm doing stars
                if stars_are_primary_type:
                    prefix = 'gas_'
                    particle_data[prefix+'SmoothingLength'] = h_all_g
                    particle_data[prefix+'Masses'] = m_all_g
                    particle_data[prefix+'Z'] = zm_all_g
                    particle_data[prefix+gasvalue_to_image] = c_all_g  
                        
                print('... ... ... sending to main imaging routine ... ... ...')               
                
                ## alright, now we can actually send this to the imaging routine 
                ##   to make the frame image! 
                image_maker.image_maker(sdir, j_of_frame, image_key=image_key, 
                  field_of_view=field_of_view, centering=np.zeros(3),
                  output_directory=output_directory, preloaded_particle_data=particle_data, 
                  saturation_level=saturation_level, **kwargs)
                
                ## don't need to do anything else here, the image & datafile are dumped!
                print('----- frame {} of {} complete! ------'.format(
                    j_of_frame, frame_number_grid[-1]), flush=True)

        print(":::::: finished with snapshot {} ::::::".format(
            snapshot_list[i_snap]), flush=True)

    return ## success!



def set_saturation_from_fov(field_of_view, rescale=1.0):
    return 0.1 * 1.e11/1.0e10 * (2.0/field_of_view)**(0.3) * rescale

def ok_scan(input,xmax=1.0e10,pos=0):
    if (pos==1):
        return (np.isnan(input)==False) & (np.abs(input)<=xmax) & (input > 0.);
    else:
        return (np.isnan(input)==False) & (np.abs(input)<=xmax);

def mask_existing_frames(frames_to_do, sdir, output_directory, 
    field_of_view, image_key, verbose=False, **kwargs):
    """
    checks if output images exist so that we can skip them.  only the
    frames that don't yet exist.
    """

    no_file_exists = np.ones(frames_to_do.size, dtype=bool)
    for kk, frame_number in enumerate(frames_to_do):
        # note that here we use the frame number as a stand-in for the snapshot number
        image_filename = image_maker.generate_directory_and_filename(sdir, 
            frame_number, output_directory=output_directory, 
            field_of_view=field_of_view, image_key=image_key, 
            **kwargs)
        if os.path.isfile(image_filename+'.pdf'):
            if verbose:
                print("skipping frame {} because it already exists in {}".format(
                    frame_number, output_directory))
            no_file_exists[kk] = False
    return frames_to_do[no_file_exists]


def idmatch(id_i, id_f): 
    ## match items in set (i) to those in set (f)
    ##  --- this is particularly speedy when the ids are all unique (case here) --
    
    # SGK -- that's not true here -- idmatch_sort is much faster and uses much less memory

    """
    In [67]: ids = np.arange(1e5)

    In [68]: ids1 = ids[np.random.permutation(ids.size)]

    In [69]: ids2 = ids[np.random.permutation(ids.size)]

    In [70]: %timeit idmatch(ids1, ids2)
    /usr/local/bin/ipython:8: DeprecationWarning: recommend using idmatch_sort instead
      if __name__ == '__main__':
    316 ms ± 37.2 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

    In [71]: %timeit idmatch_sort(ids1, ids2)
    216 ms ± 28.9 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

    and it just gets worse with bigger datasets:
    In [75]: ids = np.arange(1e6)

    In [76]: ids1 = ids[np.random.permutation(ids.size)]

    In [77]: ids2 = ids[np.random.permutation(ids.size)]

    In [78]: %timeit idmatch(ids1, ids2)
    4.39 s ± 150 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

    In [79]: %timeit idmatch_sort(ids1, ids2)
    3.5 s ± 104 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

    and they give the same answer:
    In [80]: old_return = idmatch(ids1, ids2)
    In [81]: new_return = idmatch_sort(ids1, ids2)
    In [88]: for ii in range(len(old_return)):
        ...:     print((old_return[ii] == new_return[ii]).all())
        ...:
    True
    True
    """

    import warnings
    warnings.warn("recommend using idmatch_sort instead", DeprecationWarning)

    index_dict_i = dict((k,i) for i,k in enumerate(id_i));
    index_dict_f = dict((k,i) for i,k in enumerate(id_f));
    inter = set(id_i).intersection(set(id_f));
    indices_i = np.array([ index_dict_i[x] for x in inter ]);
    indices_f = np.array([ index_dict_f[x] for x in inter ]);
    return indices_i, indices_f;

def idmatch_sort(id_i, id_f, check=False):
    """
    better way of matching indices.  does so by sorting; seems to
    be about 40% faster and use less than half the memory -- those
    dictionaries are memory hogs. 

    this works by searching for the indices that intersect in sorted
    versions of the arrays, then returning the indices that will give
    you those intersecting values (in sorted order) in the unsorted 
    arrays
    """

    import numpy as np

    #only search for these indices:
    intersect = np.intersect1d(id_i, id_f)

    #build the sorting arrays
    sorti_i = np.argsort(id_i)
    sorti_f = np.argsort(id_f)

    #gives the indices that returns the intersecting ids
    indices_i = np.searchsorted(id_i, intersect, sorter=sorti_i)
    indices_f = np.searchsorted(id_f, intersect, sorter=sorti_f)
    
    #but in order for these to be useful, I want to return indices 
    # that give you the sorted values when applied to the original array
    # so, I apply these selectors to sorti
    indices_i = sorti_i[indices_i]
    indices_f = sorti_f[indices_f]

    if check:
        assert (id_i[indices_i] == id_f[indices_f]).all()
        
    return indices_i, indices_f


def compile_matched_ids( id_i, id_f ):
    sub_id_i, sub_id_f = idmatch_sort(id_i, id_f)
    
    nomatch_from_i_in_f = (id_i > -1.0e40) ## should return all true
    if (sub_id_i.size > 0):
        nomatch_from_i_in_f[sub_id_i] = False ## is matched
    
    nomatch_from_f_in_i = (id_f > -1.0e40) 
    if (sub_id_f.size > 0):
        nomatch_from_f_in_i[sub_id_f] = False
    
    return sub_id_i, sub_id_f, nomatch_from_i_in_f, nomatch_from_f_in_i


def interpolation_for_movie_pos( dt, delta_t, u_i,u_f, vu_i,vu_f, \
    sub_id_i,sub_id_f, nomatch_from_i_in_f,nomatch_from_f_in_i , periodic=0, order=3 ):

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
        if (periodic==1):
            ## u is periodic in 2pi, so estimate its final position and wrap coordinate
            ## guess where its going to be to get 'best estimate' of the number of revolutions
            #x_exp = xi + 0.5*(vf+vi)*delta_t # actually can get some errors when vf/i switch signs...
            x_exp = xi + vi*delta_t
            ## now pin the final location to the appropriate value
            xf = get_closest_periodic_value( x_exp, xf , xi)
        
        if (order==3):
            x2 = 3.*(xf-xi) - (2.*vi+vf)*delta_t
            x3 = -2.*(xf-xi) + (vi+vf)*delta_t
            ## third-order interpolation: enables exact matching of x, v, but can 'overshoot'
            u[0:Nmatched] = xi + vi*dt + x2*tau*tau + x3*tau*tau*tau
        elif (order==2):
            ## second-order: minimizes absolute velocity difference (more stable, no exact velocity matching)
            x2 = (vf-vi)*delta_t/2.
            x1 = (xf-xi) - x2
            u[0:Nmatched] = xi + x1*tau + x2*tau*tau
        else:
            ## use linear interpolation below if for some reason this is unstable: least accurate, but most stable
            u[0:Nmatched] = u_i[sub_id_i] + (u_f[sub_id_f]-u_i[sub_id_i])*tau 
        
    ## now unmatched: first those in (i) with no match in (f)
    if(Nunmatched_i>0):
        u[Nmatched:Nmatched+Nunmatched_i] = \
            u_i[nomatch_from_i_in_f]+vu_i[nomatch_from_i_in_f]*dt;
    ## now unmatched: now those in (f) with no match in (i)
    if(Nunmatched_f>0):
        u[Nmatched+Nunmatched_i:Nmatched+Nunmatched_i+Nunmatched_f] = \
            u_f[nomatch_from_f_in_i]+vu_f[nomatch_from_f_in_i]*(dt-delta_t);

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

    return j_hat;



def interpolation_for_movie_scalar(dt, delta_t, u_i,u_f, \
    sub_id_i,sub_id_f, nomatch_from_i_in_f,nomatch_from_f_in_i, \
    mass_rescale=0, hsml_rescale=0, age_rescale=0, age_at_end=0 ):

    dt=np.array(dt); delta_t=np.array(delta_t); ## make sure casting is ok
    tau = dt/delta_t;
    Nmatched=u_i[sub_id_i].size;
    Nunmatched_i=u_i[nomatch_from_i_in_f].size;
    Nunmatched_f=u_f[nomatch_from_f_in_i].size;
    u=np.zeros(Nmatched+Nunmatched_i+Nunmatched_f);
    
    # first compute for objects in (i) with a match in (f)
    if (Nmatched>0):
        u[0:Nmatched] = u_i[sub_id_i] + (u_f[sub_id_f]-u_i[sub_id_i])*tau
    ## now unmatched: first those in (i) with no match in (f)
    if(Nunmatched_i>0):
        if(mass_rescale==1):
            u_t = u_i[nomatch_from_i_in_f] * (1.-tau**(1./5.))**5.
        elif(hsml_rescale==1):
            u_t = u_i[nomatch_from_i_in_f] / ((1.-tau**(1./3.))**3.+1.0e-5)
        elif(age_rescale==1):
            u_t = u_i[nomatch_from_i_in_f]+delta_t*tau
        else:
            u_t = u_i[nomatch_from_i_in_f]
        u[Nmatched:Nmatched+Nunmatched_i] = u_t;
    ## now unmatched: now those in (f) with no match in (i)
    if np.size(age_at_end)>1:
        age_now = age_at_end[nomatch_from_f_in_i] + dt - delta_t
        delta_age = 5.0e-4 # ...star formation timescale, in Gyr? ..ramping time..
        eta = 1./(1.+np.exp(-age_now/delta_age))
    else:
        if mass_rescale==1: eta =(1.-(1.-tau)**(1./5.))**5.
        if hsml_rescale==1: eta =(1.-(1.-tau)**(1./3.))**3.
    
    if(Nunmatched_f>0):
        if(mass_rescale==1):
            u_t = u_f[nomatch_from_f_in_i] * eta
        elif(hsml_rescale==1):
            u_t = u_f[nomatch_from_f_in_i] / (eta+1.0e-5)
        elif(age_rescale==1):
            u_t = u_f[nomatch_from_f_in_i] + dt - delta_t
            u_t[u_t<0.0]=0.0 #...try not to break it with <0
        else:
            u_t = u_f[nomatch_from_f_in_i]
        u[Nmatched+Nunmatched_i:Nmatched+Nunmatched_i+Nunmatched_f] = u_t;

    return u;


def cosmological_time(a,h=0.71,Omega_M=0.27):
    ## exact solution for a flat universe
    x=Omega_M/(1.-Omega_M) / (a*a*a);
    t=(2./(3.*np.sqrt(1.-Omega_M))) * np.log( np.sqrt(x) / (-1. + np.sqrt(1.+x)) );
    t *= 13.777 * (0.71/h); ## in Gyr
    return t;
    
    
def cross(x,y):
    #return np.cross(x,y,axis=1)
    c=0.*x
    c[:,0] = x[:,1]*y[:,2] - x[:,2]*y[:,1]
    c[:,1] =-x[:,0]*y[:,2] + x[:,2]*y[:,0]
    c[:,2] = x[:,0]*y[:,1] - x[:,1]*y[:,0]
    return c


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

    return r, v_r, phi, v_phi, j_hat

    
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
    
    return dx_corrected + x_initial


def interpolation_for_movie(dt, delta_t, \
        m_i,h_i,c_i,zm_i,\
        m_f,h_f,c_f,zm_f, \
        x_in,y_in,z_in,vx_in,vy_in,vz_in, \
        x_fn,y_fn,z_fn,vx_fn,vy_fn,vz_fn, \
        s_i,s_f,nm_i_f,nm_f_i, \
        cen_i,cen_f , use_polar_interpolation=1, stars=0 ):
    
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

    if (use_polar_interpolation==1):
        r_i,vr_i,phi_i,vphi_i,j_hat_i = xyz_to_polar(x_i,y_i,z_i,vx_i,vy_i,vz_i)
        r_f,vr_f,phi_f,vphi_f,j_hat_f = xyz_to_polar(x_f,y_f,z_f,vx_f,vy_f,vz_f)

        # r gets treated like a normal spatial coordinate (except cannot go <0; ignore for now?)
        r = interpolation_for_movie_pos(dt,delta_t, r_i,r_f,vr_i,vr_f, s_i,s_f,nm_i_f,nm_f_i,order=1)
        # predict the 'final' phi, theta, knowing that these are actually periodic:
        # note for phi, need third-order: inexact matching can leave 'wedges' missing from orbits
        phi = interpolation_for_movie_pos(dt,delta_t, phi_i,phi_f,vphi_i,vphi_f, s_i,s_f,nm_i_f,nm_f_i, periodic=1,order=3)
        x_t = r * np.cos(phi);  y_t = r * np.sin(phi);

        j_hat = interpolation_for_movie_jhat(dt,delta_t, j_hat_i,j_hat_f, s_i,s_f,nm_i_f,nm_f_i )
        ## define the absolute 'phi' relative to an arbitrary (but fixed) 90-degree rotation of the angular momentum vector j
        x_jhat = 0.*j_hat
        x_jhat[:,0]=0.*j_hat[:,0]; x_jhat[:,1]=-j_hat[:,2]; x_jhat[:,2]=j_hat[:,1]
        jj = np.sqrt( x_jhat[:,0]*x_jhat[:,0] + x_jhat[:,1]*x_jhat[:,1] + x_jhat[:,2]*x_jhat[:,2] )
        for j in [0,1,2]: x_jhat[:,j] /= jj
        ## generate y-vector by cross-product a x b = c (gaurantees right-hand rule)
        y_jhat = cross( j_hat, x_jhat )
        jj = np.sqrt( y_jhat[:,0]*y_jhat[:,0] + y_jhat[:,1]*y_jhat[:,1] + y_jhat[:,2]*y_jhat[:,2] )
        for j in [0,1,2]: y_jhat[:,j] /= jj

        x = x_t*x_jhat[:,0] + y_t*y_jhat[:,0]
        y = x_t*x_jhat[:,1] + y_t*y_jhat[:,1]
        z = x_t*x_jhat[:,2] + y_t*y_jhat[:,2]

        ## use cartesian interpolation for some points
        x_c=interpolation_for_movie_pos(dt,delta_t, x_i,x_f,vx_i,vx_f, s_i,s_f,nm_i_f,nm_f_i,order=1)
        y_c=interpolation_for_movie_pos(dt,delta_t, y_i,y_f,vy_i,vy_f, s_i,s_f,nm_i_f,nm_f_i,order=1)
        z_c=interpolation_for_movie_pos(dt,delta_t, z_i,z_f,vz_i,vz_f, s_i,s_f,nm_i_f,nm_f_i,order=1)

        ## combine velocities into initial and final matched vectors
        Nmatched=x_i[s_i].size;
        Nunmatched_i=x_i[nm_i_f].size;
        Nunmatched_f=x_f[nm_f_i].size;
        u=np.zeros(Nmatched+Nunmatched_i+Nunmatched_f);
        vr=0.*u; vphi=0.*u; rr=0.*u; rvphi_i=r_i*vphi_i; rvphi_f=r_f*vphi_f;
        if (Nmatched>0):
            vr[0:Nmatched] = np.sqrt((vr_i[s_i]**2.+vr_f[s_f]**2.)/2.)
            vphi[0:Nmatched] = np.sqrt((rvphi_i[s_i]**2.+rvphi_f[s_f]**2.)/2.)
            rr[0:Nmatched] = np.sqrt((r_i[s_i]**2.+r_f[s_f]**2.)/2.)
        if(Nunmatched_i>0):
            vr[Nmatched:Nmatched+Nunmatched_i] = vr_i[nm_i_f]
            vphi[Nmatched:Nmatched+Nunmatched_i] = rvphi_i[nm_i_f]
            rr[Nmatched:Nmatched+Nunmatched_i] = r_i[nm_i_f]
        if(Nunmatched_f>0):
            vr[Nmatched+Nunmatched_i:Nmatched+Nunmatched_i+Nunmatched_f] = vr_f[nm_f_i]
            vphi[Nmatched+Nunmatched_i:Nmatched+Nunmatched_i+Nunmatched_f] = rvphi_f[nm_f_i]
            rr[Nmatched+Nunmatched_i:Nmatched+Nunmatched_i+Nunmatched_f] = r_f[nm_f_i]
        ## attempt to guess what's 'rotation supported' and not, and use appropriate interpolation
        vr2 = vr*vr
        v2 = vr2 + vphi*vphi
        #no_rot_support = ( vr2 > v2/3. ) | ( rr > 5. ) # usually does ok, too much 'breathing' at late times
        no_rot_support = ( vr2 > v2/2. ) | ( rr > 30. ) 
        x[no_rot_support] = x_c[no_rot_support]
        y[no_rot_support] = y_c[no_rot_support]
        z[no_rot_support] = z_c[no_rot_support]

    else:
        x=interpolation_for_movie_pos(dt,delta_t, x_i,x_f,vx_i,vx_f, s_i,s_f,nm_i_f,nm_f_i)
        y=interpolation_for_movie_pos(dt,delta_t, y_i,y_f,vy_i,vy_f, s_i,s_f,nm_i_f,nm_f_i)
        z=interpolation_for_movie_pos(dt,delta_t, z_i,z_f,vz_i,vz_f, s_i,s_f,nm_i_f,nm_f_i)

    x = x + (cen_f[0]-cen_i[0]) * dt / delta_t + cen_i[0]
    y = y + (cen_f[1]-cen_i[1]) * dt / delta_t + cen_i[1]
    z = z + (cen_f[2]-cen_i[2]) * dt / delta_t + cen_i[2]
    # want to account for stellar age for star particles that fall in nomatch_f_in_i ...
    if not(stars):
        m=interpolation_for_movie_scalar(dt,delta_t, m_i,m_f, s_i,s_f,nm_i_f,nm_f_i, mass_rescale=1)
        h=interpolation_for_movie_scalar(dt,delta_t, h_i,h_f, s_i,s_f,nm_i_f,nm_f_i, hsml_rescale=1)
        c=interpolation_for_movie_scalar(dt,delta_t, c_i,c_f, s_i,s_f,nm_i_f,nm_f_i)
        #c=np.log10(interpolation_for_movie_scalar(dt,delta_t, 10.**c_i,10.**c_f, s_i,s_f,nm_i_f,nm_f_i))
        zm=interpolation_for_movie_scalar(dt,delta_t, zm_i,zm_f, s_i,s_f,nm_i_f,nm_f_i)
    else:
        # if stars.. interpolate for their formation!!
        m=interpolation_for_movie_scalar(dt,delta_t, m_i,m_f, s_i,s_f,nm_i_f,nm_f_i, mass_rescale=1, age_at_end=c_f)
        h=interpolation_for_movie_scalar(dt,delta_t, h_i,h_f, s_i,s_f,nm_i_f,nm_f_i, hsml_rescale=1, age_at_end=c_f)
        c=interpolation_for_movie_scalar(dt,delta_t, c_i,c_f, s_i,s_f,nm_i_f,nm_f_i, age_rescale=1, age_at_end=c_f)
        #c=np.log10(interpolation_for_movie_scalar(dt,delta_t, 10.**c_i,10.**c_f, s_i,s_f,nm_i_f,nm_f_i))
        zm=interpolation_for_movie_scalar(dt,delta_t, zm_i,zm_f, s_i,s_f,nm_i_f,nm_f_i, age_at_end=c_f)

    return x,y,z,m,h,c,zm;


## here's the grunt work of loading the data we'll need 
def load_pos_vel_etc(sdir, snum, GAS=0, SWAP_TEMP_RHO=0, 
    BOX_MAX=0, CENTER_BOX=[0.,0.,0.], min_stellar_age=0, 
    convert_float32=True, **kwargs):
    """
    load the particle data needed for a movie.  returns coordinates
    in physical (unless use_physical_coordinates==False is in kwargs)
    and NOT shifted to the center (but it does throw out particles 
    that are far from the center)
    """


    dum = np.zeros(1)
    
    if GAS:
        ptype = 0
    else:
        ptype = 4

    xyz_all = load_fire_snap('Coordinates', ptype, sdir, snum, **kwargs)
    if (np.array(xyz_all) == 0).all():
        # no data for the primary particle type in this snapshot
        return dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum;

    x_all, y_all, z_all = xyz_all.T
    vx_all, vy_all, vz_all = load_fire_snap('Velocities', ptype, sdir, snum, **kwargs).T

    m_all = load_fire_snap('Masses', ptype, sdir, snum, **kwargs)
    zm_all = load_fire_snap('Z', ptype, sdir, snum, **kwargs)
    h_all = load_fire_snap('SmoothingLength', ptype, sdir, snum, **kwargs)
    id_all = load_fire_snap('ParticleIDs', ptype, sdir, snum, **kwargs)

    if GAS and SWAP_TEMP_RHO:
        c_all = load_fire_snap('Density', ptype, sdir, snum, **kwargs)

        # in this case, we also set the metallicity / dust content of hot gas to zero
        temp = load_fire_snap('Temperature', ptype, sdir, snum, **kwargs)
        zm_all[temp>1e6] = 1e-6

    elif GAS:
        c_all = load_fire_snap('Temperature', ptype, sdir, snum, **kwargs)
    else:
        c_all = np.maximum(min_stellar_age, load_fire_snap('StellarAgeGyr', ptype, sdir, snum, **kwargs))

    if m_all.size <= 1: 
        return dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum;

    if convert_float32:
        myfloat = np.float32
    else:
        myfloat = np.float64

    x_all = x_all.astype(myfloat)
    y_all = y_all.astype(myfloat)
    z_all = z_all.astype(myfloat)
    
    vx_all = vx_all.astype(myfloat)
    vy_all = vy_all.astype(myfloat)
    vz_all = vz_all.astype(myfloat)

    h_all=np.array(h_all,dtype=myfloat)
    c_all=np.array(c_all,dtype=myfloat) 
    m_all=np.array(m_all,dtype=myfloat)
    zm_all=np.array(zm_all,dtype=myfloat)
    
    if (BOX_MAX != 0):
        x=x_all-CENTER_BOX[0]; y=y_all-CENTER_BOX[1]; z=z_all-CENTER_BOX[2];
        ## make sure everything is cast correctly
                
        ok=ok_scan(c_all,pos=1) & ok_scan(m_all,pos=1) & ok_scan(h_all,pos=1) & ok_scan(zm_all,pos=1) & \
            ok_scan(x,xmax=BOX_MAX) & ok_scan(y,xmax=BOX_MAX) & ok_scan(z,xmax=BOX_MAX) & \
            (h_all < 100.) & (m_all > 1.0e-12);
        
        for vec in [id_all, m_all, x_all, y_all, z_all, vx_all, vy_all, vz_all, h_all, c_all, zm_all]:
            vec = vec[ok];

    h_all *= 1.25
    return id_all, m_all, x_all, y_all, z_all, vx_all, vy_all, vz_all, h_all, c_all, zm_all


def build_snapshot_list(sdir, **kwargs):
    i=0; imax=610; snums=[-1];
    while (i < imax):
        npart_total = load_fire_snap('NumPart_Total', -1, sdir, i, no_missing_print_message=True, **kwargs)
        if (np.array(npart_total) != 0).any():
            snums.append(i)
            if (i > imax-50): #keep going if we haven't run out of files
                imax += 50;
        i += 1;
    snums = np.array(snums)
    snums=snums[snums>=0];
    return snums
    

def get_precalc_zoom_center(center_filename, time_desired, cosmological=True):
    time, x0, y0, z0 = np.loadtxt(center_filename, unpack=True)
    
    # put in comoving coords for interpolation -- will correct back to phys later
    if cosmological:
        for vec in [x0, y0, z0]:    vec /= time

    #nsmoo = 21
    #windowfun = 'hanning' ## gaussian-like smoothing (should use larger window function)

    nsmoo = 13
    windowfun = 'flat' ## boxcar smoothing

    linear_interp = False
    if time.size < 4 * nsmoo:
        print("-- not enough centers to smooth; interpolating between raw values")
        linear_interp = True

    cen = np.zeros(3)
    for j , p in enumerate([x0,y0,z0]):
        if linear_interp:
            pp = p
        else:
            pp = util.smooth(p, window_len=nsmoo, window=windowfun)

        if cosmological: 
            pp_interp = time_desired * np.interp( np.log(time_desired), np.log(time), pp )
        else:
            pp_interp = np.interp( time_desired, time, pp )
        cen[j] = pp_interp

    print("center at time = {:.2f} is ({:.2f}, {:.2f}, {:.2f})".format(time_desired, *cen))
    return cen
    

def ensure_center_of_time_exists(snapshot_list, sdir, center_filename,
    force_rebuild=False, ptype=4, **kwargs):
    """
    check if sdir + center_filename exists; if not, create it
    by looping over the snapshot_list and calling estimate_zoom_center 
    on each in turn, centering on ptype 4 (stars)
    """
    if os.path.isfile(center_filename) and not force_rebuild:
        try:
            scalefactors, x, y, z = np.loadtxt(center_filename, unpack=True)
            centers = np.vstack((x, y, z)).T
            assert scalefactors.size >= snapshot_list.size - 10  #make sure we have positions for most of the snapshots
            print("-- will use centers in {}".format(center_filename))
            return
        except Exception:
            print("-- file in {} exists, but failed to load enough centers from it; will recompute...".format(center_filename))

    # if we get here, we need to build the list
    from visualization.mm_centering import build_camera_centering
    print("-- computing centers for {} snapshots".format(snapshot_list.size))
    build_camera_centering(snapshot_list, sdir, center_filename, ptype=ptype, **kwargs)
    print("-- saved centers to {}".format(center_filename))
        



# builds the list of times to dump frames
def build_time_frame_grid( sdir, snapshot_list, frames_per_gyr,
        time_min=0, time_max=0, cosmological=0, **kwargs):

    ## set times for frames
    if(time_min==0):
        time_min = load_from_snapshot('Time', -1, sdir, snapshot_list[0], **kwargs)
    if(time_max==0):
        time_max = load_from_snapshot('Time', -1, sdir, snapshot_list[-1], **kwargs)
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

    return time_frame_grid, a_scale_grid


if __name__ == '__main__':
    """
    command line interface to movie_maker
    """

    movie_maker_defaults = util.get_default_args(movie_maker)

    import argparse
    description = "Call movie_maker on a simulation with the specified arguments."
    parser = argparse.ArgumentParser(description="Call movie_maker on a simulation.  Can add any other ")

    ### required arguments
    parser.add_argument('sdir', help="path to directory containing snapshot files")
    parser.add_argument('output_directory', help="path to directory to save the output images in")
    
    ## optional arguments
    parser.add_argument('--field_of_view', help="scale of the box size to visualize", type=float,
        default=movie_maker_defaults['field_of_view'])
    parser.add_argument('--center_filename', help="name of the file to read or save center position track, within sdir",
        default=movie_maker_defaults['center_filename'])
    parser.add_argument('--image_key', help="type of images to create", default='star')
    parser.add_argument('--frames_per_gyr', help="number of frames to make per Gyr of time passing", 
        type=float, default=movie_maker_defaults['frames_per_gyr'])

    ### arguments related to what frames/times/snapshots to visualize
    parser.add_argument('--time_min', help='initial time of the movie', type=float,
        default=movie_maker_defaults['time_min'])
    parser.add_argument('--time_max', help="final time of the movie (ignored if <=0)", type=float, 
        default=movie_maker_defaults['time_max'])

    parser.add_argument('--frame_min', type=float, default=movie_maker_defaults['frame_min'], 
        help="first frame number to process (specific to number of snaps, frames_per_gyr, etc., so be careful!)")
    parser.add_argument('--frame_max', type=float, default=movie_maker_defaults['frame_max'], 
        help="last frame number to process (ignored if <= 0)")

    parser.add_argument('--snap_min', type=float, default=movie_maker_defaults['snap_min'],
        help="first snapshot to process and include in the frame numbering.")
    parser.add_argument('--snap_max', type=float, default=movie_maker_defaults['snap_max'],
        help="last snapshot to process (ignored if <=0)")

    ### less used arguments
    parser.add_argument('--saturation_level', type=float, default=movie_maker_defaults['saturation_level'],
        help="saturation level for the movie.  if >0 (recommended, otherwise it'll change from " +
            "frame to frame), then it's in code units.  if set to 'None', 'none', or None, then it's " +
            "set according formula used in old movie maker (proportional to field_of_view**-3")

    parser.add_argument('--comoving_units', action='store_true', default=False,
        help="remove the expansion of the universe from the movie (untested!)")

    parser.add_argument('--overwrite', action='store_true', help="Overwrite existing files", 
        default=False)

    # parse the arguments (or raise an error if we didn't get what we need) into a dictionary
    args = vars(parser.parse_args())

    movie_maker(**args)



