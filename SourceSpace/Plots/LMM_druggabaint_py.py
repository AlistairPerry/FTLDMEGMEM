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


perform_LOOCV = True 


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


    #Save  
    outfname = FigOutDir + "/" + "DrugGABAInt_" + "Scatt_" + ROInames[ROI] + ".tiff"
    plt.savefig(outfname, dpi=300, format='tiff')
    plt.show()
    

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



"""

Optional ML

"""

while perform_LOOCV == True:

    """

    Perform Random Forest and LOOCV to determine performance
    
    Only right auditory cortex for computational reasons
    
    Import required libraries here

    """

    from sklearn.model_selection import LeaveOneOut
    from sklearn.model_selection import cross_val_score
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import cross_val_predict
    
    #Start
    #Load and define variables again
    dat_fname = ''.join([base_fname, 'RAUD', '.txt'])
    
    MMNtab = pd.read_table(dat_fname)
    
    
    MMNtab_vint = MMNtab[['GABA', 'Diff', 'Age']]
        
        
    #GABA and Diff values defined by 1st and 4th columns, respectively
    X = MMNtab_vint[['GABA', 'Age']]
    y = MMNtab_vint['Diff']
    
    #X = np.array(X).reshape(-1,1)
    
    #Create LOOCV procedure
    CrossVal = LeaveOneOut()
    
    #Create model
    RF_model = RandomForestRegressor(random_state=1)
    
    #Evaluate model
    CrossVal_scores = cross_val_score(RF_model, X, y, scoring='neg_mean_squared_error', cv=CrossVal)

    #Predict scores
    CrossVal_predy = cross_val_predict(RF_model, X, y, cv=CrossVal)

    corres = pearsonr(y, CrossVal_predy)

    #Force positive
    #CrossVal_scores = np.abs(CrossVal_scores)
    
    #Report performance
    print('Mean squared error is: %.3f (%.3f)' % (np.mean(CrossVal_scores), np.std(CrossVal_scores)))
    
    #Plot actual-predicted relationship
    sns.regplot(x=y, y=CrossVal_predy, color="r")
    
    #Plot attributes
    plt.xlabel(r'Actual values ($\Delta$ MEG Drug Response) ', fontsize=12, fontweight="bold")
    plt.ylabel("Predicted", fontsize=12, fontweight="bold")
    plt.title("Random Forest + LOOCV", fontsize=15, fontweight="bold")

    
    #And don't forget model results
    rval=round(corres[0], 2)
    pval=round(corres[1], 3)
    figstr = f'r={rval}, p={pval}'
    plt.text(min(y), 0.05, figstr, fontweight="bold")
    
    
    mse=round(np.mean(CrossVal_scores), 3)
    mse_sd=round(np.std(CrossVal_scores), 3)
    figstr = f'MSE={mse} (SD={mse_sd})'
    plt.text(min(y), 0.025, figstr, fontweight="bold")

    #Save  
    outfname = FigOutDir + "/" + "DrugGABAInt_" + "LOOCV_" + "RAUD" + ".png"
    plt.savefig(outfname, dpi=100, format='png')
    plt.show()
    
    #Complete
    break