#!/bin/bash

## If you plan to use this, it will require some editing, but this provides a basis for mapping your volumetric images to FreeSurfer surfaces

# Script to run surface processing.
# This script assumes you have run recon-all completely on a subject
# Then you have data stored somewhere else (DATA_DIR_IN) that you want mapped onto the surface

# Made by Christopher Rowley 2023

## Note subject 5 has no ihMT, so they have been removed from this. 

# This is the directory holding the MRI data you wish to map
DATA_DIR_IN="/data/ImageLocation"

declare -a SUBJECT=('hc01' 'hc02' 'hc03' 'hc04')

### Subject names in the FreeSurfer Subject directory "$SUBJECTS_DIR"
declare -a SUBID=('hc01' 'hc02' 'hc03' 'hc04')


## Had to redo one subject as they had differences in their prescan normalize being turned on...
#declare -a SUBJECT=('hc01')

### Subject names in the FreeSurfer Subject directory "$SUBJECTS_DIR"
#declare -a SUBID=('Nelson_HC01')

## INPUT IMAGE SUFFIX
INsuffix="mnc.gz"


CONTRASTS=(sparseMP2RAGE_M0,sparseMP2RAGE_T1,spMP2RAGE_MTsat_2k)
##### Redo just for the 2k    
#CONTRASTS=(spMP2RAGE_MTsat_2k)

#####################################################################################################################################
## END INPUTS!!! LET THIS RUN
#####################################################################################################################################
ITER=0
for i in "${SUBJECT[@]}"
do
    # Export to each subjects native directory
    DATA_DIR_OUT=$DATA_DIR_IN/$i/minc/freeSurfer/ants
    #IMGDIR=$DATA_DIR_IN/$i/minc/register/matlab
    IMGDIR=$DATA_DIR_IN/$i/minc/ants/matlab
    
    # Put the ants results in a new folder to allow me to compare after. 
    mkdir $DATA_DIR_OUT
    
    for j in ${CONTRASTS//,/ }
    do
        #echo $j
        #echo "$DATA_DIR_OUT"
        #echo "$IMGDIR"
        echo "${SUBID[$ITER]}"
        printf "\n"
        

        # First need to convert to nifti
        mnc2nii $IMGDIR/$j.$INsuffix $IMGDIR/$j.nii
        
        ## Sampling mid depth use: --projfrac 0.5
        ## Average across cortical depths: --projfrac-avg 0.2 0.8 0.2

        # Leave out smoothing, can do that in matlab, if you want to do it here, add flag ( --surf-fwhm 3  )
        mri_vol2surf --src $IMGDIR/$j.nii --out $DATA_DIR_OUT/$j.lh.mgh --regheader ${SUBID[$ITER]} --hemi lh --surf white --projfrac-avg 0.4 0.6 0.2 --trgsubject fsaverage
        mri_vol2surf --src $IMGDIR/$j.nii --out $DATA_DIR_OUT/$j.rh.mgh --regheader ${SUBID[$ITER]} --hemi rh --surf white --projfrac-avg 0.4 0.6 0.2 --trgsubject fsaverage

        # Remove the nifti file for space saving
        rm $IMGDIR/$j.nii

    done
    
    IMGDIR=$DATA_DIR_IN/$i/minc/register/matlab
    
    ## Then do the B1 map, but leave separate since its a difference folder
    # First need to convert to nifti
    mnc2nii $IMGDIR/../resampled_b1field.mnc $IMGDIR/b1field.nii

    # Leave out smoothing, can do that in matlab, if you want to do it here, add flag ( --surf-fwhm 3  )
    mri_vol2surf --src $IMGDIR/b1field.nii --out $DATA_DIR_OUT/b1field.lh.mgh --regheader ${SUBID[$ITER]} --hemi lh --surf white --projfrac-avg 0.4 0.6 0.2 --trgsubject fsaverage
    mri_vol2surf --src $IMGDIR/b1field.nii --out $DATA_DIR_OUT/b1field.rh.mgh --regheader ${SUBID[$ITER]} --hemi rh --surf white --projfrac-avg 0.4 0.6 0.2 --trgsubject fsaverage

    # Remove the nifti file for space saving
    rm $IMGDIR/b1field.nii
    
    
    ITER=$(expr $ITER + 1)
    #echo $ITER

done


        
# Here with this flag you can take the average over the cortex in steps of 0.2
#--projfrac-avg 0 1 0.2













