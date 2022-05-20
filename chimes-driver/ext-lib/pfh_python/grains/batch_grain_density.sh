#!/bin/sh

procmin=0
procmax=16

# Obtain our unique rank: 0 - (numprocs-1)
proc=$1   # get from parent job array
#proc=${PBS_VNODENUM}

# these keys matter only if we start switching the mode
SNAPDIR="'/scratch/01799/phopkins/grain_turbulence/rho0pt01mach10/output'"

n0=0;n1=0;
case $proc in
0) n0=0;n1=5; ;;
1) n0=6;n1=10; ;;
2) n0=11;n1=15; ;;
3) n0=16;n1=20; ;;
4) n0=21;n1=25; ;;
5) n0=26;n1=30; ;;
6) n0=31;n1=35; ;;
7) n0=36;n1=40; ;;
8) n0=41;n1=45; ;;
9) n0=46;n1=50; ;;
10) n0=51;n1=55; ;;
11) n0=56;n1=60; ;;
12) n0=61;n1=65; ;;
esac

runf=grain_density.${proc}.run
echo import grains.grain_density_from_snapshot as gd > ${runf}
echo x=gd.grain_density_from_snapshot_loop\(sdir=$SNAPDIR,snum_min=$n0,snum_max=$n1\) >> ${runf}

if (( $proc >= $procmin && $proc <= $procmax)); then
python ${runf} > ${runf}.g.log
fi
