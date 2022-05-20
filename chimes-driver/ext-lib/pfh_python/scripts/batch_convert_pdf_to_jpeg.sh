#!/bin/sh

for i in $1*
do
  echo $i 
  ## split off the filename without the extension 
  pathfilename=${i%.*} 
  ## extension 
  fileext=${i##*.} 
  convert $pathfilename.pdf -quality 100 -resize 800x800 $pathfilename.jpeg 
done