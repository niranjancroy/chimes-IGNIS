#!/usr/bin/env python3

#SBATCH --job-name=m12i-star-movie
#SBATCH --output=movie/star.log

### please change these!
#SBATCH --mail-user=sheagk@gmail.com
#SBATCH --mail-type=all

### Note:  for a typical m12X_res7100 run (e.g. m12i_res7100),
### the memory usage of a typical process is <~ 15 GB.  Keep 
### this number in mind when setting ntasks-per-node below.
### For example, on the Rusty nodes, I can 

#### typical Stampede options:
#SBATCH --partition=skx-normal
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=8
#SBATCH --time=48:00:00
#SBATCH --account=TG-AST140064

#### typical stampede development/testing options
##SBATCH --partition=skx-dev
##SBATCH --nodes=2
##SBATCH --ntasks-per-node=8
##SBATCH --time=2:00:00
##SBATCH --account=TG-AST140064

### typical Rusty options:
##SBATCH --partition=cca
##SBATCH --nodes=4
##SBATCH --ntasks-per-node=14
##SBATCH --exclusive
#### if you need even more memory, ensure that you're using the Skylake
#### nodes, which have 728 GB RAM and 40 cores.  Otherwise, you'll probably
#### end up using Broadwell node(s), which have 512 GB Ram and 28 cores.
###SBATCH --constraint=skylake

from subprocess import call
import numpy as np
import os
import sys
import shutil
from glob import glob

import visualization.movie_maker as movie_maker
import pfh_utils as util
from gizmopy.load_fire_snap import evaluate_if_cosmological

################################################## SET VARIABLES ################################################
# snapshot directory -- i.e. where snapshot_xxx.hdf5 or snapdir_xxx/ live:
sdir = 'output'

# output directory -- i.e. where the images and log files are saved
output_directory = 'movies/stars_40kpc'

# field of view
field_of_view = 40

# name of the file to load or create with the center position 
# as a function of scale factor (in physical units)
center_filename = 'movie_center.txt'

# what type of movie to make?
image_key = 'star'

## any other options to pass, e.g. frames_per_gyr, whatever, as '--key=val' items
other_options = []    

## set one of these to point to either:
## a custom installed launcher (set the directory)
hardcoded_launcher_dir = None
## or a custom installed disBatch (set the disBatch.py path)
hardcoded_disbatch_command = None


################## these don't need to by reset typically ##################
## what's our python command?
python = 'python3'

## where do we find movie_maker.py?
script = os.environ['PYDIR']+'/visualization/movie_maker.py'
################################################ end of setting vars #############################################

## make sure the output directory exsits
os.makedirs(output_directory, exist_ok=True)

## turn the skip list into a flat array, if it isn't already
# frame_skip_list = np.array(frame_skip_list).flatten()

## make sure that calls were set correctly (if at all):
if hardcoded_launcher_dir is not None:
    assert os.path.isdir(hardcoded_launcher_dir), "you set 'launcher' directory by hand, but it doesn't exist!"
if hardcoded_disbatch_command is not None:
    assert os.path.isfile(hardcoded_disbatch_command), "you set 'disBatch.py' command by hand, but it doesn't exist!"
    assert hardcoded_disbatch_command.endswith('disBatch.py'), "hardcoded_disbatch_command should end with 'disBatch.py'"

######## figure out what our minimum and maximum frame numbers are:
## check if frames_per_gyr, snap_min/_max, and/or time_min/_max are set in other_options
movie_maker_defaults = util.get_default_args(movie_maker.movie_maker)
def get_ops_or_fallback_to_default(key):
    ops_list = list(filter(lambda x: key in x, other_options))
    val = movie_maker_defaults[key] if len(ops_list) == 0 else float(ops_list[0].split('=')[-1])
    return val

frames_per_gyr = get_ops_or_fallback_to_default('frames_per_gyr')
snap_min = get_ops_or_fallback_to_default('snap_min')
snap_max = get_ops_or_fallback_to_default('snap_max')
time_min = get_ops_or_fallback_to_default('time_min')
time_max = get_ops_or_fallback_to_default('time_max')

### copying the code from movie maker that builds the frame number grid...
snapshot_list = movie_maker.build_snapshot_list(sdir)
if snap_max <= 0: 
    snap_max = np.max(snapshot_list)

time_frame_grid, a_scale_grid = movie_maker.build_time_frame_grid(sdir, 
    snapshot_list, frames_per_gyr, time_min=time_min, time_max=time_max, 
    cosmological=evaluate_if_cosmological(sdir, snapshot_list[0]))

## and now I have what I need to figure out my expected frame numbers:
frame_number_grid = np.arange(time_frame_grid.size)
n_orig = frame_number_grid.size

## now figure out which ones I already have:
frame_number_grid = movie_maker.mask_existing_frames(frame_number_grid, 
    sdir, output_directory, field_of_view, image_key)
n_todo = frame_number_grid.size
print("skipping {} frames because they already exist".format(
    n_orig - n_todo), flush=True)

def split_list_of_numbers_with_skip(numbers, toskip, num_blocks, warning_size=2):
    """
    split up a list of integers going into num_blocks contiguous blocks,
    skipping the integers listed in toskip.

    tries to perform the most even split possible, but obviously can
    be rendered tricky by the items in toskip

    warns if any of the blocks are <= warning_size (but still continues)

    returns:
        two arrays of min/max values, as strings
    """

    numbers = np.array(numbers, dtype=int)
    numbers_todo = np.setdiff1d(numbers, np.array(toskip, dtype=int))

    ### now have to find contiguous blocks:
    # don't want to do this with the itertools+operator recipe in utils because
    # that ignores the fact that I have a maximum block size set by wanting to 
    # share the work among all the available tasks
    max_block_size = (numbers_todo.size / num_blocks) + 1  #just about even split
    blocks = [] 

    start_idx = 0
    idx = 1
    while idx < numbers_todo.size:
        if (numbers_todo[idx] - numbers_todo[start_idx] != idx - start_idx) or idx - start_idx > max_block_size:
            # then I'm no longer in a continguous block or I've built a block too big
            # go back a step and build a block
            end_idx = idx - 1

            ends = [ numbers_todo[start_idx], numbers_todo[end_idx] ]
            if ends[1] - ends[0] <= 0:
                print("-- warning: found a block that's {} to {}".format(ends[0], ends[1]))
            blocks.append(ends)
            start_idx = end_idx  #make the blocks overlap by one

        #eitehr way, increment the index
        idx += 1

    #now deal with the final block:
    end_idx = idx - 1
    ends = [ numbers_todo[start_idx], numbers_todo[end_idx] ]
    if ends[1] - ends[0] <= warning_size:
        print("-- warning: found a block that's {} to {}".format(ends[0], ends[1]))
    blocks.append(ends)

    nblocks = len(blocks)
    block_sizes = [b[1]-b[0] for b in blocks]
    print("Split into {} blocks with median sizes {} (range={} to {})".format(nblocks, np.median(block_sizes), min(block_sizes), max(block_sizes)))

    #now pick out the longest num_blocks blocks, if there are too many:
    if nblocks > num_blocks:
        print("Subselecting the {} largest blocks.".format(num_blocks))
        sorti = np.argsort(block_sizes)[::-1]   #go from largest to smallest
        blocks = np.array(blocks)              
        my_blocks = blocks[sorti[:num_blocks]]     #pick out the n largest
        print("-- warning:  skipped the following blocks because there were too many splits:")
        for block in blocks:
            if block not in my_blocks:
                print('{} - {}'.format(*block))
        blocks = my_blocks

    block_beginnings = np.array([str(b[0]) for b in blocks])
    block_ends = np.array([str(b[1]) for b in blocks])
    return block_beginnings, block_ends

def get_logfilename(basename):
    assert type(basename) == str
    if not os.path.isfile(basename):
        return basename
    counter = 1
    while True:
        if not os.path.isfile(basename+'.'+str(counter)):
            return basename + '.' + str(counter)
        counter += 1


ntasks = int(os.environ.get('SLURM_NTASKS', 1))
jobid = os.environ.get('SLURM_JOB_ID', 'DUMMY')

##### now split the snapshots to do and build the launcher jobfile:
jobfname = 'movie_commands-'+jobid+'.sh'
commands = []

## split the list of frames to do into min-max pairs
## don't need to worry about those that we're skipping because they were already removed
frame_min_ar, frame_max_ar = split_list_of_numbers_with_skip(
    frame_number_grid, [], ntasks)

## one command per min-max pair:
for jj in range(frame_min_ar.size):
    fmin = frame_min_ar[jj]
    fmax = frame_max_ar[jj]
    logfn = get_logfilename(output_directory + '/log.' + image_key + 
        '-'+ str(field_of_view) + 'kpc-' + fmin + '-' + fmax) 

    # create as a list, then smash together at the end
    cmd = [ python, script, sdir, output_directory, 
        '--center_filename='+center_filename, '--image_key='+image_key, 
        '--field_of_view={}'.format(field_of_view), 
        '--frame_min={}'.format(fmin), '--frame_max={}'.format(fmax)]
    cmd += other_options

    ## join together with a space and append our log filename to the end
    cmd = ' '.join(cmd) + ' &>'+logfn
    commands.append(cmd)

num_jobs = len(commands)

#write all that to a file that we can call with either launcher, disbatch, or parallel:
with open(jobfname, 'w') as f:
    f.write('\n'.join(commands)+'\n')

print("Wrote {} commands to {}".format(len(commands), jobfname))

#### now call one of our three options for running the commands:
call_launcher = False
call_disbatch = False


## first check if we set either of our options by hand:
if hardcoded_launcher_dir is not None:
    call_launcher = True
    launchdir = hardcoded_launcher_dir
elif hardcoded_disbatch_command is not None:
    call_disbatch = True
    disbatch_command = hardcoded_disbatch_command
elif 'LAUNCHER_DIR' in os.environ:
    call_launcher = True
    launchdir = os.environ['LAUNCHER_DIR']
elif shutil.which('disBatch.py') is not None:
    call_disbatch = True
    disbatch_command = 'disBatch.py'  # should be in path if we loaded module
else:
    ## TODO fall back to GNU parallel here...
    print("Can't load either launcher or disBatch")
    print("Please load the appropriate module or install it and set the path in this file")
    sys.exit(1)

### figure out how to split these up across the nodes:
num_nodes = int(os.environ['SLURM_JOB_NUM_NODES'])
cpus_per_node = int(os.environ['SLURM_CPUS_ON_NODE'])

## OK, so I have a  total of num_jobs split over num_nodes => num_jobs / num_nodes jobs per node
ntasks_per_node = num_jobs/num_nodes  #don't round this; round later

## but I have cpus_per_node per node, so that means that I have (cpus/node) /(jobs/node) cpus per job
cpus_per_task = int(np.floor(cpus_per_node / ntasks_per_node))   #round this down

## now round ntasks_per_node
ntasks_per_node = int(np.floor(ntasks_per_node))

## this won't do anything if not using the openmp version fof the raytrace
## also won't do anything if openmp thread count set by hand
## but if not, it'll automatically try to use all the threads
os.environ['OMP_NUM_THREADS'] = str(cpus_per_task)


## now actually make the call:
if call_launcher:
    ## call launcher to split up the work:
    os.environ['LAUNCHER_DIR'] = launchdir
    os.environ['LAUNCHER_PLUGIN_DIR'] = launchdir+'/plugins'

    ## launcher reads all the cpus / task etc from environment variables
    ## we just have to make sure it knows that we're using SLURM
    os.environ['LAUNCHER_RMI'] = 'SLURM'

    os.environ['LAUNCHER_WORKDIR'] = os.getcwd()
    os.environ['LAUNCHER_JOB_FILE'] = jobfname

    print("Calling launcher from {}".format(launchdir))
    call([launchdir+'/paramrun'])

elif call_disbatch:
    ## call disBatch to split up the work
    print("calling {}".format(disbatch_command))

    dbatch_logdir = output_directory + '/disbatch-logs/'
    os.makedirs(dbatch_logdir, exist_ok=True)
    call([disbatch_command, 
            '-p', dbatch_logdir,   ## store log files etc. in the output directory
        jobfname])

else:
    print("don't know how we got here, but no way to call the command file...")
    sys.exit(1)