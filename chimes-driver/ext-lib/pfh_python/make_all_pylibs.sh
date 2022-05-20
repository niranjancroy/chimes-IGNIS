#!/bin/sh

cwd=$(pwd)
#if [[ $cwd = *"pfh_python"* ]]; then
if [ -z $PYDIR ]; then
    PYDIR=$cwd
    echo "!!!! --- warning ---- !!!!"
    echo "!!!! --- warning ---- !!!!"
    echo "Setting PYDIR to the present dir; $PYDIR"
    echo "Add the following line to your .bashrc or lots of things won't work:"
    echo "export PYDIR=$PYDIR"
fi

PYLIBDIRS=(
$PYDIR/agn_spectrum
$PYDIR/attenuation
$PYDIR/c_libraries/LOS_column_singlePOV
$PYDIR/c_libraries/RayTrace_RGB
$PYDIR/c_libraries/RayTrace_RGB_Serial
$PYDIR/c_libraries/SmoothedProjPFH
$PYDIR/c_libraries/StellarHsml
$PYDIR/c_libraries/hop_halofinder
$PYDIR/c_libraries/fof
$PYDIR/c_libraries/neighbor_finder
$PYDIR/c_libraries/CorrFn
)
for i in ${PYLIBDIRS[*]}
do
  echo 'making...'$i
  cd $i
  make clean
  make
done
