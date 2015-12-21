#!/bin/sh

\rm -r log
for j in 01 02 03  
do
  echo "Processing file " $j
  mkdir log
  ../../src/cac 2 image.pgm mask_${j}.pgm cage_${j}.txt
  cd log
  ../../genimages 640 480 | gnuplot 
  cd ..
  mv log log_${j}
done

