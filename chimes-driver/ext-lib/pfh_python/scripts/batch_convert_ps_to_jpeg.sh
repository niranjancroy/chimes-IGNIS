#!/bin/sh

for i in $1*
do
  echo $i
  
  ## split off the filename without the extension
	pathfilename=${i%.*}
  ## extension
	fileext=${i##*.}
  
  
  #pstogif $i 
  #convert $i -quality 100 -resize 800x800 frame0350.jp2
  
  convert $pathfilename.ps -quality 100 -resize 800x800 $pathfilename.jpeg

done

