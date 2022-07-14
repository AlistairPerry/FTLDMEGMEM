# -*- coding: utf-8 -*-
"""
LMM_drug_gaba_int

This script file replicates the R function in this repository

For examining the influence of GABA on patient MEG responses to memantine

"""


"""
Setup

"""

#Import paths

import pandas as pd

import os 


#Set data dir

work_dir = "/Users/alistairperry/Documents/Cambridge/Project/Analysis/LFPs/newmaxfilter_icafixes/LFPs_COH_wrad_allMEGch_wfids/C_Plots/DrugAnalysis/MRS"

os.chdir(work_dir)


#Regions to load from text file

base_fname = 'LFPs_MRSIFGgabacorassoc_meanMMN3_DrugDiff_Pats_fullsample_stats_LMMtableforR'

ROInames = ['RIFG', 'RSTG', 'RAUD']


"""
Start for loop, picking up each ROI in every loop

"""

for ROI in ROInames:
    
    dat_fname = base_fname + ROInames[ROI] + '_.txt'

    print(dat_fname)