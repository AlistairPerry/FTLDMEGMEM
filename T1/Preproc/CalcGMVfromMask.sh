#!/bin/bash

DATADIR=$1

OUTdir=${DATADIR}/GMV_4849OP8

mkdir ${OUTdir}


mask=/imaging/rowe/Michelle/AlistairIFG/A_Scripts/rs4849OP8.nii


#Note for full IFG area = #voxels 2933
 
find ${DATADIR}/ -name "mwc1m*" > normfiles.txt


while read line
do


subjfname=$(basename $line)


#Remove norm prefix and extract subject

basefname=$(echo $subjfname | cut -c5- ) 

subj=$(echo $basefname | cut -d "_" -f 1)


subj_lc=$(echo "$subj" | awk '{print tolower($0)}')



#Now lets get to calc the GMV

#fslmeants -i $line -m ${mask} -o ${OUTdir}/${subj_lc}_meanGM_fullIFG.txt

fslstats -K ${mask} ${line} -V -M | awk {'print $1 * $3'} > ${OUTdir}/${subj_lc}_GMV_fullIFG.txt


#Sure we can do calc mean in same script - nah this is fine for now


done < normfiles.txt
