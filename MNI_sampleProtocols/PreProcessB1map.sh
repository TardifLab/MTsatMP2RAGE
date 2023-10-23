#!/bin/bash

# Sample code to process thee B1 map, and resample it.
# This assumes that your files are in minc format... If not, you
# should register the anatomical B1 image to the T1w image (rigid)
# then apply that transform to the phase image.

#######################################################################3

cd /image/directory
REF=anat-mpm_acq-mni_megre_T1w_1mm.mnc.gz
LOWRES=fmap-b1_tfl1.mnc.gz
    
mincresample -clobber -float -like $REF -fill $LOWRES temp.mnc 

minccalc -float -expr "clamp(A[0]/800,0,3)" temp.mnc B1map_rs.mnc

rm temp.mnc
    
