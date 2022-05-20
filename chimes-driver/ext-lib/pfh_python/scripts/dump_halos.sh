#!/bin/sh
#proc=$1
proc=${PBS_VNODENUM}

SDIR0="'/scratch/m/murray/phopkins/"
sdir1="m12_mr_tst'"
sdir2="m12_mr'"
sdir3="zoom_nosn'"
minmass=1.0e7
case $proc in
1) SDIR=${sdir1}; ns=441; ;;
2) SDIR=${sdir2}; ns=441; ;;
3) SDIR=${sdir3}; ns=441; ;;
4) SDIR=${sdir1}; ns=190; ;;
5) SDIR=${sdir1}; ns=240; ;;
6) SDIR=${sdir1}; ns=290; ;;
7) SDIR=${sdir1}; ns=340; ;;
8) SDIR=${sdir1}; ns=440; ;;
esac

SDIR2=${SDIR0}${SDIR}
runf=dump_halos_${proc}.py
echo import halo_finder > ${runf}
echo halo_finder.compile_halo_properties\($SDIR2,$ns,min_halomass_res=$minmass/1.e10\) >> ${runf}
python ${runf} > ${runf}.log

