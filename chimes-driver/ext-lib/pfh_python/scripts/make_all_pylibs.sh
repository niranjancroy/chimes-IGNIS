#!/bin/sh

PYLIBDIRS=(
$PYDIR/agn_spectrum
$PYDIR/attenuation
$PYDIR/c_libraries/LOS_column_singlePOV
$PYDIR/c_libraries/RayTrace_RGB
$PYDIR/c_libraries/SmoothedProjPFH
$PYDIR/c_libraries/StellarHsml
$PYDIR/c_libraries/hop_halofinder
$PYDIR/c_libraries/fof
)
for i in ${PYLIBDIRS[*]}
do
  echo 'making...'$i
  cd $i
  make clean
  make
done
