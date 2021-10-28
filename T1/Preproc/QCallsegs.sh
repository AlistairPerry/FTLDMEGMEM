#!/bin/bash
#Loop through all subjects seg images to assess quality


DATADIR=$1

cd ${DATADIR}


#Do separately for controls and patients

#Controls


find ${DATADIR}/* -name "C*denoised.nii" > /imaging/rowe/Michelle/AlistairIFG/A_Scripts/ControlSubjs.txt


while read line

do

T1fname=$(basename $line)

freeview -v $T1fname c1${T1fname}:heatscale=0.1,1:colormap=Heat:opacity=0.2 -layout 1 -viewport y -zoom 1.25


done < /imaging/rowe/Michelle/AlistairIFG/A_Scripts/ControlSubjs.txt



#Patients


find ${DATADIR}/* -name "P*denoised.nii" > /imaging/rowe/Michelle/AlistairIFG/A_Scripts/PatientSubjs.txt


while read line

do

T1fname=$(basename $line)

freeview -v $T1fname c1${T1fname}:heatscale=0.1,1:colormap=Heat:opacity=0.2 -layout 1 -viewport y -zoom 1.25


done < /imaging/rowe/Michelle/AlistairIFG/A_Scripts/PatientSubjs.txt


#Finito
