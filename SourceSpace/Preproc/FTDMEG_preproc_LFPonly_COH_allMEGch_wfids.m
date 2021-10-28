function FTDMEG_preproc_LFPonly_COH_allMEGch_wfids(subj, subjfname, subjt1scan, subjDIR)

%Coregister individuals MRI & MEG or use template MRI if no MRI available 
%(here, 'S' is the structural MRI location for this individual):
%PP_Coreg_wfids([subjDIR '/' subjfname], subj, subjt1scan)
PP_Coreg_wfids([subjDIR '/' subjfname], subj, subjt1scan)

%Create an inversion state for the subsequent extraction of sources from the
%sensor space data. This is currently set up for roving mismatch negativity
%trials, and uses only MEGPLANAR sensors. If you want to alter this, you can
%do so directly by altering the script:
PP_Inversion_COH_allMEGch([subjDIR '/' subjfname])

%Create directory so LFPs can be sourced easily
mkdir([subjDIR '/' 'LFPs_COH_allMEGch_wfids'])

%Extract the 6 sources of interest to your directory of choice:
PP_Extract([subjDIR '/' subjfname], [subjDIR '/' 'LFPs_COH_allMEGch_wfids/'])

end