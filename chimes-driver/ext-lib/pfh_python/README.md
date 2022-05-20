# Visualization instructions:

While this repo contains a great many tools, the majority of users will be interested in either the `gizmopy` sub-package or the `visualization` sub-package.  Instructions for the former are the main gizmo repository (where the code is also available); I describe how to use the latter here.

Before getting started with visualizations, please, please, please read and follow all of the following instructions.  We've tried to reduce the complexity of the calls, but there are still several moving parts involved in creating an image/movie out of these tools, and it can take some effort to get a nice-looking image out.  Moreover, don't be afraid to turn to Photoshop/GIMP at a certain point -- there's no shame in editing the levels of an image by hand to make it look good (though obviously this is more difficult for movies).

If you used the previous version of these tools, then feel free to [skip to the section below](#quickstart-guide-for-users-of-the-previous-version) where I summarize what's changed in this new version.

## Setup:
### Requirements:
* Python 3.  This may work with python 2, but I doubt it.
* `scipy`, `numpy`, `h5py`, and `matplotlib`       
* An OpenMP-enabled C compiler.  Stampede2's Intel compiler will work fine.  For now, the compilation step checks for `icc`, `gcc-8`, and `gcc-7` before defaulting to the system `gcc`.
* This repository (`pfh_python`)
* One of the following to spread the `movie_maker` calls over many CPUs/nodes (not needed if just making images or planning to spread calls yourself):
    * the TACC `launcher` utility.  `launcher` is accessible with `module load launcher` on Stampede2, but should be installable in your user directory on other computers too (though I haven't tried).  Installation instructions are available at [the Github repository](https://github.com/TACC/launcher).  Note that you'll have to edit the `submit_pfh_movie.py` script discussed below to point to your personal installation.

    * the Flatiron `disBatch` utility.  `disBatch` is accessible with `module load disBatch` on Rusty, and there are again instructions for installing on other computers [in the Github repository](https://github.com/flatironinstitute/disBatch).    Note that you'll again have to edit the `submit_pfh_movie.py` script discussed below to point to your personal installation.

### Setup of your environment
Make sure that `python3` works and allows you to import `h5py`, `scipy`, `matplotlib`, and `numpy`.  The first can be a problem on Stampede2 where their default build seems to be broken -- I've had success creating a virtual environment and installing `h5py` (and the other reqs) in that environment, rather than using the system build, but I also understand that you should be able to get everything working with just a local install of these modules (i.e., after loading the `python3` module, do `python3 -m pip install --upgrade --user numpy scipy h5py matplotlib`)

#### Setup this repository:
> **WARNING:** the `setup.py` file included in this repository is **NOT** a `Distutils` `setup.py` file. That is, `python setup.py` will **NOT** compile and install the code.  Please follow the instructions below. 

This repository includes Python scripts that frequently call underlying C libraries.  In order to use those C libraries, you must tell Python where to find them and compile them into executables.  

The former is accomplished via the `PYDIR` environment variable, which should point to the root directory of this repository.  In other words, you should be able to do `ls $PYDIR/visualization/movie_maker.py` and find the `movie_maker.py` file.  This must be set (correctly) both when you compile the code and whenever you want to run it, so I recommend adding it to your `.bashrc`.  You'll also need `PYDIR` in your `PYTHONPATH` in order to let Python know where to import the  scripts from.  The compilation script will tell how to do this if you run it without first setting `PYDIR`.

So, we'll clone the repository, set `PYDIR`, add `PYDIR` to our `PYTHONPATH`, then compile all the C libraries using the included script that runs `make` on each underlying folder.  I'm going to assume that this repository is going to live in `~/code`, and that you've already loaded the appropriate module for your compiler.  Nonetheless, these commands may need to be modified slightly depending on your environment:

        mkdir -p ~/code/
        cd ~/code
        
        # clone the repository:
        hg clone https://bitbucket.org/phopkins/pfh_python/
        
        # set PYDIR and add it to your PYTHONPATH
        echo "export PYDIR=$HOME/code/pfh_python" >> ~/.bashrc
        echo "export PYTHONPATH=$PYDIR:$PYTHONPATH" >> ~/.bashrc
        source ~/.bashrc
        
        # cd into the directory and run the script to run make on all the libraries
        cd $PYDIR
        ./make_all_pylibs.sh

On the first run of the final command, there'll be some errors due to the script trying to remove old versions of the libraries that don't exist yet, but you should mostly see a bunch of `gcc` commands, with an `icc` command mixed in with the appropriate openmp flags.  To be safe, I often rerun `make_all_pylibs.sh`, which should eliminate all of the errors.

You should now be able to start up `python` and do, e.g., `from visualization.image_maker import image_maker`.  If you can't, check your `PYTHONPATH`.

*** 

## Making an image:
The work-horse of making a movie is making a bunch of images, so it's important to make sure that you can successfully do that.  The image making is done by the function `image_maker` in `visualization.image_maker.py`.  That function will (by default) load up a snapshot, compute stellar smoothing lengths based on the distances to some number of nearest neighbor particles (if making an image of the stars), do the ray-trace in three bands, turn those maps into an image, and save the image and the raw maps.

In order to save time on subsequent runs, the code will attempt to save the stellar softening lengths in an HDF5 file that, by default, goes in the snapshot directory.  It therefore helps to have write access to that directory.  If you don't, the code will still run, but it'll have to compute those smoothing lengths every time you want to make an image of that snapshot.

So, let's try making a still image of a single snapshot.  I **STRONGLY** recommend reading the docstring of `image_maker.image_maker`, but the most important arguments (in my opinion) are: 

* `sdir` (required), the snapshot directory.  This should contain either snapshot_xxx.hdf5 or snapdir_xxx/snapshot_xxx.y.hdf5.  For example, if you're in the simulation directory, it should probably be 'output'.
* `snum` (required), the snapshot number you want to visualize (e.g. 600)
* `image_key`, which sets the type of image that you'll be making.  Set this to 'star' to create the standard mock-color star-light images and 'gas_Temperature' to create the standard gas images.  See the docstring of `visualization.image_maker.image_maker` for more details and options (e.g. you can set `image_key='star_Temperature'` to get a stellar image with a 1-color gas temperature map overlaid.)
* `dynamic_range`, which sets the dynamic range of the image.  I recommend setting this somewhat high, since you can always bring it down in post with GIMP.
* `field_of_view`, which sets the physical size of the image.
* `pixels`, which sets the resolution of the image.

Putting all that together, we might run the following in an `ipython` session:

    from visualization.image_maker import image_maker
    image = image_maker('output',   # setting sdir = 'output'
                        600,        # setting snum = 600 -- either 'output/snapshot_600.hdf5' or 'output/snapdir_600' *must* exist 
                        image_key='star', # make an image of the starlight
                        image_name_prefix='test_image',  #beginning of output image name (has some other data appended on)
                        field_of_view=50., # 50 kpc-across image
                        pixels=2048)

Once it finishes, it should tell you that it saved a something.  Add `.pdf` to that filename to see the image that you created (and make sure to zoom in -- though it may be small, it's high resolution and can be resized up.)

You might also want to play around with the following arguments to `image_maker` (w/ default values listed):
    
* `saturation_level=-0.0001` Sets either the absolute saturation level (if positive) or the fraction of pixels that should saturated (if negative).  e.g. -0.0001 says that 0.01% of the pixels will be saturated.
* `centering=''` How to center the image.  if a blank string, will use the largest object in the box.  set to a length 3 array to set the center by hand.  should be in physical code units.
* `projection=''` What sort of projection to do.  if 'camera' is not in the string, then you get a flat field.  if projection=='camera', then a camera is placed at the center position and pointed along the z-axis (with pitch/roll/yaw controlled via the euler_angles arguement).  if projection is either 'camera_tocenter' or 'camera_to_center', then the camera will be shifted along the z-axis so that it points *at* the center position.  see the docstring for more info.
* `camera_opening_angle=90`  The opening angle for the camera, if there is one.  larger values give more a fish-eye effect.
* `output_directory='./'`  Where to save the output images (gets put before `image_name_prefix`)
* `labels='scale_redshift'` What labels to add to the image (i.e. whether or not to include the time/redshift and a scale bar)
* `h_rescale_factor=3` Arbitrary multiplicative scaling for the stellar smoothing lengths (doesn't impact gas at all).  Set to a larger value to spread out the light more -- makes it look less point-sourcey, but also can give strange behavior at early times when the universe is rapidly expanding

In general, the "look" of the image is basically determined by the maximum and dynamic range of the colorbar, which are controlled by `saturation_level` and `dyanmic_range`, respectively.  Setting the saturation level dyanmically (i.e. make it less than 0) is good for single images, but doesn't really work with a movie.  If you're setting `saturation_level` manually (i.e. `saturation_level > 0`) and the image is completely black (white), then you need to lower (raise) your saturation level.  `dynamic_range` sets the dynamic range/stretch of the image -- to make less luminous pixels show up, you can raise `dynamic_range`, such that the minimum luminosity that isn't black is smaller.  Think of `dynamic_range` as setting the lowest luminosity that will be visible (i.e. not black), given that you've already fixed the maximum of the color/luminosity scale.  Of course, the correlary of a higher `dynamic_range` is that higher luminosities will be scrunched together.  

Once you run this command, you should see the routine go and load up the data (in serial, of course, but only loading the necessary properties), find the correct center, then compute the smoothing lengths for the stars (unfortunately still in serial, but this file will at least be saved and re-used), then start the ray-trace (in parallel).

The raytrace is parallelized by having each thread compute an image, then layering the images from each thread together.  That means that increasing the number of pixels or increasing the number of threads will require more RAM -- a 2048x2048 image requires about 512 MB using 8 threads to store the image.  This is pretty negligible compared to the size of the particle data, but if you set `openmp=50`, then it can get to be a problem.  In other words, I don't recommend setting `openmp` too high.  I've had good speedups without a bit memory hit with `openmp` somewhere in the range of 8 - 16.

If that worked, awesome!  You can play around with the `openmp` option to see how much of a speedup you get (remember if you set it equal to False or 0 then the code should load up the serial library).  We're now ready to make a movie, but first, a note about a very useful function for...

***

### Making edge-on/face-on images
If you're interested in making a still image of the galaxy, there's a good chance that you want to view the galaxy either *edge-on* or *face-on*.  Fortunately, there's a pre-written function that'll make this easy for you, `image_maker.edgeon_faceon_projection`, which will calculate the proper projection matrix for the view you want:
    
    from visualization.image_maker import edgeon_faceon_projection
    edgeon_faceon_projection('outout',   #snapshot director
                             600,        #snapshot number
                             field_of_view=40,  #set the size of the image
                             faceon=True,       #do a faceon image (or set faceon=False or edgeon=True to get an edge-on image)
                             **kwargs)  #kwargs are passed to image_maker, get_center, and load_fire_snap

Note that this function can accept any valid `image_maker` keyword arguments (e.g. `projection`, `dynamic_range`, etc.)

## Making a movie:
Now that you can make an image, you should be able to make a movie!  Before we get into that though, some words about what's going to happen and how we're going to parallelize it:

In short, the routine will:

1.  loop over all the snapshots
2.  interpolate the particle properties (position, mass, age, metallicity, temperature, etc.) in between the snapshots
3.  create images using `image_maker.image_maker` on the interpolated particle properties.  

For the stars, this means loading up the stars and the gas:  the stars add luminosity to the pixel grid according to some smoothing length and the gas attenuates the luminosity behind them.  For gas movies, only the gas is required -- typically, the pixels are then colored by their temperature.  That means that the stellar movies are slightly more memory intensive -- keep this in mind when trying to figure out how to best parallelize the tasks.

Now, on the topic of parallelization, the movie making basically consists of these steps:
1.  Figure out the center in each snapshot to create a smooth track for the camera
2.  Load up snapshot 000 and 001, interpolate the particle properties between them, and make the images at each requested time that falls at/in between those snapshots.
3.  Repeat step 2 for snapshots 001 and 002.

You can probably see already that step 2 can be nearly trivially parallelized -- in practice, you could do 599 independent tasks that each go between two snapshots.  In principle, the memory usage of each task would be such that you usually can't stuff more then about 10 tasks on a single node (and sometimes even fewer), so running 599 tasks would require approximately 60 nodes -- definitely possible, but maybe not ideal.  Since there's some start-up cost before the snapshot loop, we'd also be wasting some initialization time.  Even if we do push to that many nodes, we'd be wasting an enormous fraction of the CPUs if we don't also parallelize the actual image creation.  Futhermore, the number of frames between two snapshots is not always equal, and the time it takes to trace a frame can vary as well.

I'm not accounting for the final problem, but I account for most of the rest I lay out above by parallelizing the ray-tracing with OpenMP (see above) and make it easy to parallize the loop over the snapshots/frames.  The latter uses the `launcher` utility if it's available then looks for `disBatch` if it's not.  In principle, it should be able to fall back to GNU `parallel`, but I don't have a lot of experience working with it, so it's not working right now.

***

So, without further ado, let's dive into making a movie:

#### Computing the centers:

The `movie_maker` routine requires that the center positions be computed before starting the loop over snapshots in order to create a smooth camera path.  It looks for this path in a 4 column text file located at `center_filename` (one of the arguments to `movie_maker`).  If no file exists at that path, then `movie_maker` will create it using the default technique (see below).  If that method doesn't give good results, you can try one of the other three functions that are built in, or you can define it yourself.

> If you create the center path by hand, then it should be a 4 column text file that gives (in order) the scale factor (or, for non-cosmological runs, the time), then the x-y-z position of the center _in physical code units_ (for FIRE, this should be physical kpc; i.e. not physical kpc/h, comoving kpc, or comoving kpc/h).

There are four functions that can create the center in `visualization.mm_centering`.  Two leverage Andrew Wetzel's routines, and the other two use built-in techniques.  Each option has a method for going forward in time (i.e. using center at t_{i-1} as a guess) and going backwards in time (i.e. using the center at t_{i+1} as a guess).  They are:

* `build_camera_centering`:  this is the default centering algorithm.  It starts at t = 0 then goes forward in time.  At each step, it uses the previous two guesses to estimate where the center should be and how much it is moving from step-to-step, which in theory gives a speedup.  Also has support simply use the center-of-mass of a given particle type, to use the center of mass of all black hole particles, or to use a single fixed center at all times.  

* `build_camera_centering_backwards`:  this is a variant of the default algorithm that starts at the final timestep and works backwards from there.  Most of the options are the same as `build_camera_centering`, though it lacks support for the simpler centering methods (for which there's no reason to work backwards).  However, it does have an option to pass in a `final_position_guess`, which starts off the backwards chain
    on an object near that position.

* `build_camera_centering_wetzel`:  this function leverages the center finding algorithms in Andrew Wetzel's `utilities` package, and uses the `gizmo_analysis` package to quickly read the data in the appropriate format.  It moves forward in time from the first snapshot, and uses the center from the previous time as a guess for the current time (but doesn't allow for an initial guess at t = 0).  It can also handle multiple hosts (though it only does one host at a time) by setting `host_index > 0` (i.e. to find the path of the second host in an LG run, set `host_index = 1`).  My tests suggest that this method is relatively prone to failure if you start too early, so I don't recommend using it.

* `build_camera_centering_wetzel_backwards`:  Just like `build_camera_centering_wetzel`, this leverages `utilities` and `gizmo_analysis` to do the center finding, but here it starts at the final timestep and works backwards from there.  That means that you can set `final_position_guess` to the position of the host(s) in the final snapshot.  Note, however, that `final_position_guess` must be a list of centers equal to the number of hosts that you're looking for (i.e. if `host_index==1`, then `final_position_guess` should be either None (to let the code figure it out itself), or a len-2 list of size-3 arrays).

To reiterate, either:
    
1. Skip the centering step and let `movie_maker` do it for you using `build_camera_centering`
2. Call one of the four functions above by hand ahead of time
3. Create a 4 column text file some other way

To use either `build_camera_centering` or `build_camera_centering_backwards`, do something like:

    from visualization.mm_centering import build_camera_centering
    snapshot_numbers = np.arange(601)  ### if doing all snapshots
    build_camera_centering(snapshot_numbers, 'output', 'movie_center.txt')

or:

    from visualization.mm_centering import build_camera_centering_backwards
    snapshot_numbers = np.arange(601)  ### if doing all snapshots
    build_camera_centering_backwards(snapshot_numbers, 'output', 'movie_center.txt', 
        final_position_guess=position_in_final_snapshot)

To use either of the `wetzel` routines, make sure that `utilities` and `gizmo_analysis` are installed or in your PYTHONPATH, then do something like:

    from visualization.mm_centering import build_camera_centering_wetzel
    snapshot_numbers = np.arange(601)
    build_camera_centering_wetzel(snapshot_numbers, 'output', 'movie_center.txt')

or

    from visualization.mm_centering import build_camera_centering_wetzel_backwards
    snapshot_numbers = np.arange(601)
    build_camera_centering_wetzel(snapshot_numbers, 'output', 'movie_center.txt', 
        final_position_guess=position_in_final_snapshot or None)


#### Submitting a parallel loop over the snapshots

Once you have a track for the center of the target galaxy (or have decided that you're going to let `movie_maker` create it for you, you'll want to submit a job to actually create the frames.  In order to most equally share the work (since most of the work is in the image creation, not the loading and interpolation), it's best to split the work according to the frame numbers that each job should do.

I've created a script that will parallelize the frame creation by setting up a list of `movie_maker` calls, each with different values for `frame_min` and `frame_max`, then will use either `launcher` or `disBatch` (whichever is available) to run all of those jobs in parallel.  The base script is `$PYDIR/tools/submit_pfh_movie.py`, so make a copy of it to your working directory and start editing (or set up the calls yourself if you're a masochist).

You'll need to edit the following:

  1.  Change the SLURM options (i.e. number of nodes, type of nodes, and CPUs/node) to appropriate values.  For reference, I'm typically able to stuff 3 - 4 stellar movies (of a res7100 run) on a single `skx` node.  If some error out due to memory constraints, it won't mess up the rest of the jobs, so you can be somewhat aggressive.  I also recommend starting with a development queue to make sure everything is working ok.

  2. Set the variables at the top of the script to the movie that you want to create:    
      * `sdir`:  The directory with the snapshot files or snapshot directories (usually 'output')
      * `output_directory`:  This is the directory you want to save the images to.
      * `field_of_view`:    Field of view to use for the movie, typically in physical kpc (unless you use other non-standard options)
      * `center_filename`: this should be the path to the text file that contains the position as a function of scale factor, or path to where you want to save that file if it doesn't exist yet.
      * `image_key`:  type of image to make (i.e. 'star', 'gas_Temperature', etc.)
      * `other_options`: list of any other optional arguments to pass to `movie_maker` (e.g. can put '--overwrite' in this list to overwrite existing files.).  Set as either '--key=val' items or '--key' if a store_true option.  You can (but don't need to as the code will default to only using the snapshot it finds) set `time_min`, `time_max`, `snap_min`, and `snap_max` here.  However, do _not_ set `frame_min` or `frame_max` here, as the script is going to split the frames to do over the tasks that you're assigning it.  Note that you could also set '--verwrite=True' to have the code re-make existing images.  Note that the check is done before the particle data is loaded, so it's reasonable effecient at skipping completed snapshots.

  3.  If you've installed your own `launcher` or `disBatch`, set the appropriate `hardcoded_launcher_dir` (should point to the directory) or `hardcoded_disbatch_command` (should point to disBatch.py).

Now you should be good to go -- just make sure that you're environment is set up correctly and submit.  In practice, that means making sure that `which python3` returns the python3 that lets you import everything you need, then run `sbatch submit_pfh_movie.py`.

Now you get to sit back and watch the log files for a print of either the smoothing length computation or a progress bar, which'll first show up once it starts doing the ray-tracing.

#### Compiling the frames into a movie

Once the runs are all finished, you'll be left with a whole bunch of image files that you need to compile into a movie.  I recommend using `avconv`, which part of the `libav` package (installable via Homebrew).  You'll probably also need to convert the PDFs to PNGs; for that I recommend  `convert`, part of `imagemagick`.  The syntax etc. can be tricky for these, so I recommend using the included `tools/compile_frames_to_movie.py`:

    python3 $PYDIR/tools/compile_frames_to_movie.py

This script whill take care of converting the PDFs to a good size and setting the appropriate syntax for `avconv`.  The most important arguments are `-r`, which sets the frame rate (frames per second) and `-s`, which sets the size of the movie (e.g. `-s 2000x2000`).

### Other useful scripts and functions in tools/
* `identify_missing_frames.py`: As mentioned above, run this on the base of the output frames to get a list of which ones are missing.  It'll also create the list of continuous frames that are completed, which you can skip on a continuation.  However, the `submit_pfh_movie.py` script should handle this itself now, as it checks for which frames already exist.


## Quickstart guide for users of the previous version:
First, let me just list the big changes that were made in this update:

Big changes:

1.  numexpr is not (at this time) supported, so drop all use_numexpr arguments 
2.  `gizmo_analysis` and `utilities` are no longer required (though they won't cause any conflicts either).  Particle loading is now handled by the included `gizmopy` sub-module.  As a consequence, the `snapshot_times.txt` file is no longer required (though it doesn't hurt to have it!)
3.  Support for creating gas and star images together no longer exists.  I hope to add it back in at some point, but it's not a huge priority.  In practice, this means we spend roughly twice as much time as necessary (assuming you want to make both movies) on the loading and interpolating steps.
4.  Support for making the edge-on and face-on images in a single call has been removed.  Again, this effectively means roughly twice as much time spent on loading, centering, and calculating the projection vectors.
5.  No support currently exists for any sort of yt-based images or movies. Again, hoping to add that in soon, but it's not a huge priority.
6.  `test_edgeon_faceon` has been renamed to `edgeon_faceon_projection`.
7. `submit_pfh_movie.py` has been revamped.  It's simpler and only handles a single movie at a time now, but it's much cleaner in how it does that.  It also has support for either launcher or disBatch now.  Instead of putting in snapshots to skip, you should now put in the frames that are completed.

## TODO items:
* Add support in `submit_pfh_movie.py` for GNU parallel
* Change all imports to be `from pfh_python ...` to clean up namespaces.  Ideally, even make this repo into a package that can be installed with a `setup.py` (note that the existing `setup.py` is NOT that type of `setup.py`)
* Add the required hooks to make movies of the DM density field w/ `yt`
* Or, alternatively, figure out a way to deposit the DM particles onto an image, probably somewhat like the gas density.  Should be relatively trivial to do a simple density map then hand that to matplotlib to get a colorbar.

### Completed TODOs:
* Set up `submit_pfh_movie.py` to go into the output directory, figure out which frames need to be done, and then split up the frames to do instead of requiring the list of frames to skip to be put in by hand.
* Allow the user to pass in a final position and trace the object backwards in time.
