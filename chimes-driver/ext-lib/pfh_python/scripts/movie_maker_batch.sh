#!/bin/sh

## script to parallelize batch processing of movie frames: call from a master 
##    PBS/SLURM/etc script, for some number of processors, this farms out a set of 
##    frames to each processor individually

# set range of processor list below that a given call can access (if need to limit)
procmin=4
procmax=$(($procmin+2))

# Obtain our unique rank: 0 - (numprocs-1)
#proc=$1   # get from parent job array
proc=${PBS_VNODENUM} 
# (this will be master-script specific, but should correspond to 'which processor' we're on)

# these actually define the settings of the movie script, of course can add more
SNAPDIR="'$WORK/m12_mr'"  # directory containing snapshots
OUT_DIR="'$WORK/images/'" # master directory for dumping images
IMTYPE="'gas'" # flag for star/gas/xr image
BOX=50. # box size (in physical kpc)
COSMO=1 # flag for cosmological simulation
FRAME_PER_GYR=100. # frames per Gyr of movie
EXT="'_r50'"

s0=0 # minimum snapshot number (0 = all)
s1=0 # maximum snapshot number (0 = all)

# list which defines first (n0) and last (n1) frame to be farmed out to each process
case $proc in
0) n0=706;n1=782; ;;
1) n0=783;n1=860; ;;
2) n0=1012;n1=1093; ;;
3) n0=1094;n1=1175; ;;
4) n0=1176;n1=1266; ;;
5) n0=1267;n1=1357; ;;
6) n0=1050;n1=1174; ;;
7) n0=1175;n1=1299; ;;
8) n0=1300;n1=1500; ;;
esac

cd $PBS_O_WORKDIR 
# (on some systems this is needed to be sure we're in the right directory, on others not)

# ok now use the above definitions to make a run script
runf=movie_maker.${proc}.run
echo import visualization.movie_maker as mm > ${runf}
echo q=mm.movie_maker\(xmax_of_box=$BOX,snapdir=$SNAPDIR,outputdir_master=$OUT_DIR,\
add_extension=$EXT,frame_min=$n0,frame_max=$n1,i_snap_min=$s0,i_snap_max=$s1,\
show_gasstarxray=$IMTYPE,cosmological=$COSMO,frames_per_gyr=$FRAME_PER_GYR\) >> ${runf}

# and (if in the allowed process range) run that script in python
if (( $proc >= $procmin && $proc <= $procmax)); then
python ${runf} > ${runf}.g.log
fi
