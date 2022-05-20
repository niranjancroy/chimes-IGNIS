import numpy as np
import h5py as h5py
import os.path
import math
import gc
import sys


def calculate_com_center(snapdir, snapnum, ptype=4, **kwargs):
    """
    calculate the center of mass of a given particle type
    in a single snapshot

    kwargs are passed to load_fire_snap
    """
    from gizmopy.load_fire_snap import load_fire_snap
    positions = load_fire_snap('Coordinates', ptype, snapdir, snapnum, **kwargs)
    masses = load_fire_snap('Masses', ptype, snapdir, snapnum, **kwargs)

    masses /= np.sum(masses)
    return np.sum(np.multiply(positions.T, masses).T, axis=0)


def build_camera_centering_wetzel_backwards(snums, snapdir, 
    center_filename, host_index=0, species='star', 
    final_position_guess=None, **kwargs):
    """
    use Andrew Wetzel's centering algorithm to compute the center 
    for the snapnums and save to a file.  goes backwards in time, so
    you have the option to specify a guess for the center at the 
    final time.  note that you must specify a guess for all the hosts
    in the run (i.e. final_position_guess.size should be 3 * (host_index+1))

    must have gizmo_analysis and utilities installed or in your 
    PYTHONPATH to use this function

    valid kwargs are 'snapshot_name_base' (used to define ReadClass)
    and any valid arguments to ut.particle.get_center_positions other
    than part, species, center_number, and center_positions
    """
    from copy import deepcopy

    from gizmo_analysis.gizmo_io import ReadClass
    import utilities as ut

    snapshot_name_base = kwargs.pop('snapshot_name_base', 'snap*[!txt]')
    Read = ReadClass(snapshot_name_base=snapshot_name_base, verbose=False)

    ## search for the number of hosts we're indexing to
    center_number = host_index + 1

    decreasing_times = np.zeros(snums.size)
    decreasing_centers_in_comoving = np.zeros((snums.size, 3))

    ## sort the snapshot numbers to go from large to small
    decreasing_snums = np.sort(snums)[::-1]

    for ii, snapshot_number in enumerate(decreasing_snums):
        part = Read.read_snapshots(species=species, snapshot_value_kind='index',
            snapshot_values=snapshot_number, properties=['position', 'mass'], 
            snapshot_directory=snapdir+'/', assign_host_coordinates=False)

        if species not in part:
            print("-- no {} particles in snapshot {}".format(species, snapshot_number))
            continue

        scale = part.snapshot['scalefactor']
        decreasing_times[ii] = scale

        center_guess_list = None
        if ii == 0 and final_position_guess is not None:
            center_guess_list = final_position_guess
        elif ii > 0:
            ## use the guess from the previous iteration (the next time) as a center
            center_guess_list = deepcopy(this_center_position)

        this_center_position = ut.particle.get_center_positions(
            part, species=species, center_number=center_number, 
            center_positions=center_guess_list, return_array=False, 
            **kwargs)

        decreasing_centers_in_comoving[ii, :] = this_center_position[host_index]
    
    ## now swap the order, save, and return:
    times = decreasing_times[::-1]
    centers_in_comoving = decreasing_centers_in_comoving[::-1]
    centers_in_phys = (centers_in_comoving.T * times).T

    ## skip any steps where we didn't find a center because no particles (i.e. where all three are zero)
    msk = (centers_in_phys[:,0] != 0) | (centers_in_phys[:,1] != 0) | (centers_in_phys[:,2] != 0)
    np.savetxt(center_filename, np.vstack((times[msk], centers_in_phys[msk].T)).T)
    print("saved result to {}".format(center_filename))

    return times, centers_in_phys


def build_camera_centering_wetzel(snums, snapdir, center_filename, 
    host_index=0, species='star', **kwargs):
    """
    use Andrew Wetzel's centering algorithm to compute the center 
    for the snapnums and save to a file.  goes forward in time, so
    no guess is ever needed.

    must have gizmo_analysis and utilities installed or in your 
    PYTHONPATH to use this function

    valid kwargs are 'snapshot_name_base' (used to define ReadClass)
    and any valid arguments to ut.particle.get_center_positions other
    than part, species, center_number, and center_positions
    """
    from copy import deepcopy

    from gizmo_analysis.gizmo_io import ReadClass
    import utilities as ut

    snapshot_name_base = kwargs.pop('snapshot_name_base', 'snap*[!txt]')
    Read = ReadClass(snapshot_name_base=snapshot_name_base, verbose=False)

    times = np.zeros(snums.size)
    centers_in_comoving = np.zeros((snums.size, 3))

    ## search for the number of hosts we're indexing to
    center_number = host_index + 1

    ## make sure the snapshot numbers are sorted
    snums = np.sort(snums)

    this_center_position = None

    with open(center_filename, 'a') as outfi_cen:
        for sidx, snapshot_number in enumerate(snums):
            part = Read.read_snapshots(species=species, snapshot_value_kind='index',
                snapshot_values=snapshot_number, properties=['position', 'mass'], 
                snapshot_directory=snapdir+'/', assign_host_coordinates=False)

            if species not in part:
                print("-- no {} particles in snapshot {}".format(species, snapshot_number))
                continue

            scale = part.snapshot['scalefactor']
            times[sidx] = scale

            center_guess_list = None
            if this_center_position is not None:
                ## use the guess from the previous time as a center (if we have one)
                center_guess_list = deepcopy(this_center_position)

            this_center_position = ut.particle.get_center_positions(
                part, species=species, center_number=center_number, 
                center_positions=center_guess_list, return_array=False, 
                **kwargs)

            centers_in_comoving[sidx, :] = this_center_position[host_index]
            this_host_cen_in_phys = centers_in_comoving[sidx, :] * scale

            ## now write out this scale factor...
            v = np.array([ scale, this_host_cen_in_phys[0], this_host_cen_in_phys[1], this_host_cen_in_phys[2] ]).reshape(-1,4)
            np.savetxt(outfi_cen,v)
    
    ## now return the times and centers
    centers_in_phys = (centers_in_comoving.T * times).T
    return times, centers_in_phys

    


def build_camera_centering(snums, snapdir, center_filename,
    center_on_bh = False, # center on first BH particle (good for nuclear sims)
    center_on_com = False, # center on center of mass (good for merger sims)
    set_fixed_center = False, # option to set a fixed center -- pass the center
    four_char=0, skip_bh=0,
    min_searchsize=1.,search_deviation_tolerance=0.03,
    verbose=False, ptype=4, **kwargs):

    from gizmopy.load_fire_snap import evaluate_if_cosmological
    cosmological = evaluate_if_cosmological(snapdir, snums[-1], **kwargs)

    ## ok, if we're here, center-table doesn't exist, so let's build it:
    n_snums = np.array(snums).size
    time_all = np.zeros( (n_snums) )
    cen = np.zeros( (n_snums,3) )
    n_prev = 0

    center_on_zero = False
    set_fixed_center = np.asarray(set_fixed_center)
    if (set_fixed_center.size == 1) and set_fixed_center.all():
        ## if we're told to do a fixed center, but aren't given
        ## a specific center, assume we're just centering on zero
        print("Setting all centers to zero...")
        center_on_zero = True
    elif set_fixed_center.size == 3:
        print("Setting all centers to ({:.3f}, {:.3f}, {:.3f}".format(*set_fixed_center))
    
    with open(center_filename, 'a') as outfi_cen:
        for i in range(n_snums):
            print('centering for snapshot {}..'.format(snums[i]), end='')
            
            fname, fname_base, fname_ext = check_if_filename_exists(snapdir,snums[i],four_char=four_char)
            file = h5py.File(fname,'r')
            time_all[i] = file["Header"].attrs["Time"]
            hubble = file["Header"].attrs["HubbleParam"]
            file.close()
            cen[i,:] = [0.,0.,0.]

            if center_on_zero: 
                pass
            elif set_fixed_center.size==3:
                cen[i,:] = set_fixed_center
            elif center_on_com:
                ## center on the center of mass of the default particle type
                cen[i, :] = calculate_com_center(snapdir, snums[i], ptype=ptype, **kwargs)
            elif center_on_bh:
                ## center on the center of mass of the black holes (ptype = 5)
                cen[i, :] = calculate_com_center(snapdir, snums[i], ptype=5, **kwargs)
        
            else:
                ## if not using any of the methods above, use our fancier iterative 
                ##   centering solution for the 'dominant' structure
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

                cen[i,:], n_prev = fast_zoom_center( snapdir, snums[i], \
                    max_searchsize=max_searchsize, min_searchsize=min_searchsize, search_deviation_tolerance=search_deviation_tolerance,
                    cen_guess=cen_guess,ptype=ptype,four_char=four_char,cosmological=cosmological, verbose=verbose)
            
            print('..zoom_cen at ',cen[i,:], flush=True)
            
            ## all done, write out these results if it's not zero (i.e. if we found a center)
            ## or if we're told to center on zero explicitely 
            if not (cen[i] == 0).all() or center_on_zero:
                v = np.array([ time_all[i], cen[i,0], cen[i,1], cen[i,2] ]).reshape(-1,4)
                np.savetxt(outfi_cen,v)
    return time_all, cen


def build_camera_centering_backwards(snums, snapdir, center_filename,
    final_position_guess=None, four_char=0, skip_bh=0, min_searchsize=1., 
    search_deviation_tolerance=0.03, verbose=False, ptype=4, **kwargs):

    from gizmopy.load_fire_snap import evaluate_if_cosmological
    cosmological = evaluate_if_cosmological(snapdir, snums[-1], **kwargs)

    ## ok, if we're here, center-table doesn't exist, so let's build it:
    n_snums = np.array(snums).size
    time_all = np.zeros( (n_snums) )
    cen = np.zeros( (n_snums,3) )
    n_prev = 0

    decreasing_snums = np.sort(snums)[::-1]
    decreasing_centers = cen[::-1]
    decreasing_time_all = time_all[::-1]

    for i in range(n_snums):
        print('centering for snapshot {}...'.format(decreasing_snums[i]), end='')
        
        fname, fname_base, fname_ext = check_if_filename_exists(
            snapdir, decreasing_snums[i], four_char=four_char)

        with h5py.File(fname,'r') as file:
            decreasing_time_all[i] = file["Header"].attrs["Time"]
            hubble = file["Header"].attrs["HubbleParam"]

        max_searchsize = 1000.
        rcorr = 0.*decreasing_time_all + hubble
        if(cosmological==1): 
            rcorr /= decreasing_time_all

        # start with a guess of zeros (no guess), but swap to previous guess if we have one            
        cen_guess=np.zeros(3)
        if i == 0 and final_position_guess is not None:
            cen_guess = np.array(final_position_guess)

        elif i > 0:
            cen_guess = decreasing_centers[i-1,:] * rcorr[i-1]
            if i > 1 and n_prev>100.:
                ## get a sense of how much the halo has moved over the previous two timesteps
                d_cen = np.sqrt(np.sum(
                    (decreasing_centers[i-1,:]*rcorr[i-1] - 
                    decreasing_centers[i-2,:]*rcorr[i-2])**2
                    ))
                #thold = 8.*d_cen/rcorr[i]
                thold = 5.*d_cen/rcorr[i]
                if thold < max_searchsize: 
                    max_searchsize=thold
                #if(max_searchsize < 15.*min_searchsize  ): max_searchsize=15.*min_searchsize
                if max_searchsize < 10.*min_searchsize: 
                    max_searchsize=10.*min_searchsize

        ## if not using any of the methods above, use our fancier iterative 
        ##   centering solution for the 'dominant' structure
        decreasing_centers[i,:], n_prev = fast_zoom_center( 
            snapdir, decreasing_snums[i],
            max_searchsize=max_searchsize, min_searchsize=min_searchsize, 
            search_deviation_tolerance=search_deviation_tolerance,
            cen_guess=cen_guess, ptype=ptype, four_char=four_char, 
            cosmological=cosmological, verbose=verbose)
        print('zoom_cen at ',decreasing_centers[i,:], flush=True)
            
    ## all done; write out results where not zero (i.e. where a center is found)
    time_all = decreasing_time_all[::-1]
    cen = decreasing_centers[::-1]

    msk = (cen[:,0] != 0) | (cen[:,1] != 0) | (cen[:,2] != 0)
    np.savetxt(center_filename, np.vstack((time_all[msk], *cen[msk].T)).T)
    print("saved result to {}".format(center_filename))

    return time_all, cen




def fast_zoom_center(sdir, snum, snapshot_name='snapshot', extension='.hdf5', four_char=0,
    max_searchsize=1000., min_searchsize=1., search_deviation_tolerance=0.05,
    cen_guess=[0.,0.,0.], ptype=4, cosmological=1, verbose=False):
    
    fname,fname_base,fname_ext = check_if_filename_exists(sdir,snum,\
        snapshot_name=snapshot_name,extension=extension,four_char=four_char)
    cen = np.array(cen_guess)
    if(fname=='NULL'): return cen
    if(fname_ext!='.hdf5'): return cen

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
                    file.close()
            xyz_prev = xyz
        else:
            d = np.amax(np.abs(xyz_prev-cen),axis=1)
            ok = np.where(d < d_thold)
            n0 = ok[0].shape[0]
            xyz = xyz_prev.take(ok[0],axis=0)
            if(n0 > 1):
                xyz_m += np.sum(xyz,axis=0)
                n_t += n0
        niter += 1;
        if(n_t[0] <= 0): 
            d_thold *= 1.5
        else:
            xyz_m /= n_t[0]
            d_cen = np.sqrt(np.sum((cen-xyz_m)**2))
            if verbose:
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





def test():
    for snum in [32,33]:
        if(snum==32): cen=np.array([ 3455.8801012 ,  3614.45947891,  3654.45549142])
        if(snum==33): cen=np.array([ 3455.8801012 ,  3614.45947891,  3654.45549142])
        if(snum==32): cen_guess=np.array([ 3457.03729451,  3612.84589234,  3657.79841436])
        if(snum==33): cen_guess=np.array([ 3546.72016253,  3710.18996351,  3759.31516493])
        ke_zoom_center('.',snum,ptype=0,cen_guess=cen,max_searchsize=10.)

def ke_zoom_center(sdir, snum, snapshot_name='snapshot', extension='.hdf5', four_char=0,
    max_searchsize=1000., min_searchsize=1., search_deviation_tolerance=0.05,
    cen_guess=[0.,0.,0.], ptype=0, cosmological=1):
    
    fname,fname_base,fname_ext = check_if_filename_exists(sdir,snum,\
        snapshot_name=snapshot_name,extension=extension,four_char=four_char)
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
    cen = np.array(cen_guess) * rcorr
    max_searchsize *= rcorr; min_searchsize *= rcorr; search_deviation_tolerance *= rcorr;
    if(max_searchsize > boxsize): max_searchsize=boxsize
    d_thold = max_searchsize
    if(cen[0]==0.): d_thold=1.e10
    gc.collect()
    niter = 0
    xyz = np.zeros((0,3)); vxyz=np.zeros((0,3)); m=np.zeros((0))
    xyz_m=np.zeros(3); n_t=np.zeros(1); 
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
                    m = np.concatenate([m,np.array(file['PartType'+str(ptype)+'/Masses/']).take(ok[0])])
                    vxyz = np.concatenate([vxyz,np.array(file['PartType'+str(ptype)+'/Velocities/']).take(ok[0],axis=0)])
                    xyz_m += np.sum(xyz_all,axis=0)
                    n_t += n0
                gc.collect()
            file.close()
        gc.collect()
    gc.collect()
    wt = m / np.sum(m)
    vxyz *= np.sqrt(time)
    m /= hubble
    v0 = np.median(vxyz,axis=0)
    vxyz -= v0
    v2 = np.sum(vxyz*vxyz,axis=1)
    print('KE == ',np.sum(m*v2),np.sum(m),np.max(v2))



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
