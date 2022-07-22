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
import numpy as np
from scipy.stats import pearsonr
import seaborn as sns
import matplotlib.pyplot as plt
import os 


#Set data dir

work_dir = "/Users/alistairperry/Documents/Cambridge/Project/Analysis/LFPs/newmaxfilter_icafixes/LFPs_COH_wrad_allMEGch_wfids/C_Plots/DrugAnalysis/MRS"

os.chdir(work_dir)


#Regions to load from text file

base_fname = 'LFPs_MRSIFGgabacorassoc_meanMMN3_DrugDiff_Pats_fullsample_stats_LMMtableforR_'

ROInames = ['RIFG', 'RSTG', 'RAUD']


#Output
#Set figure output

FigOutDir = "/Users/alistairperry/Documents/Cambridge/Project/Analysis/LFPs/newmaxfilter_icafixes/LFPs_COH_wrad_allMEGch_wfids/C_Plots/DrugAnalysis/MRS/Scatter_py"

#Check it exists

if os.path.exists(FigOutDir):
    pass

else:
    
    os.mkdir(FigOutDir)


"""
Start for loop, picking up each ROI in every loop

"""

#Start output variables

r2_all = np.empty(len(ROInames), dtype = float)
p_all = np.empty(len(ROInames), dtype = float)


for ROI in range(len(ROInames)):
    
    #Identify filename based on parts
    
    dat_fname = ''.join([base_fname, ROInames[ROI], '.txt'])
    
    #Read txt file with pandas
    
    MMNtab = pd.read_table(dat_fname)
    
    #Perform corr
    
    corres = pearsonr(MMNtab["GABA"], MMNtab["Diff"])

    #Extract values
    
    r2_all[ROI] = corres[0]**2
    
    p_all[ROI] = corres[1]
    
    
    #Figures
    sns.regplot(data=MMNtab, x="GABA", y="Diff", color="r")    
    
    #Plot attributes
    plt.xlabel("GABA rIFG (cor)", fontsize=12, fontweight="bold")
    plt.ylabel("MMN (PLA-MEM)", fontsize=12, fontweight="bold")
    plt.title(ROInames[ROI], fontsize=15, fontweight="bold")
    plt.show()


    #Bonus: LOOCV
   

"""

Write text output

"""

#Combine output in data frame

data = {"ROI": ROInames, "r2": r2_all, "p": p_all}

df = pd.DataFrame(data)

print("The basic statistical output is: \n")

print(df)

#Writing out text

outfname = FigOutDir + "/" + "DrugGABAInt_corrvalues.csv"

df.to_csv(outfname, sep=",", index=False)