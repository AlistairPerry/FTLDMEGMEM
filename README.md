# FTLD, MEG, and MEMANTINE

This repository contains processing, analysis, and visualisation scripts that were performed in the following paper:

* Perry, A., Hughes, L., Adams, N., Naessens, M., Murley, A., Rouse, M., ... & Rowe, J. (2022). The neurophysiological effect of NMDA-R antagonism of frontotemporal lobar degeneration is conditional on individual GABA concentration. Translational Psychiatry, available at https://www.nature.com/articles/s41398-022-02114-6


<br />

This study investigated the influence of the pharmacological agent memantine on frontotemporal brain networks in persons with frontotemporal lobar degeneration (FTLD).

<br />

The resources contained here can be allocated into 3 sections:

1. [Processing pipeline for MEG data](https://github.com/AlistairPerry/FTLDMEGMEM/tree/main/SourceSpace/Preproc)
2. [Analysis calculation of MEG responses, and group and drug differences](https://github.com/AlistairPerry/FTLDMEGMEM/tree/main/SourceSpace/Analysis)
3. [Code to reproduce publication plots in R, Python, and Markdown](https://github.com/AlistairPerry/FTLDMEGMEM/tree/main/SourceSpace/Plots)

<br />
<br />


Most importantly, the following information details how the main findings were produced, relating to sections 2/3.

<br />

## R

These scripts will produce the following findings:

1. [MEG responses across placebo and mematine sessions for control and patient populations](https://github.com/AlistairPerry/FTLDMEGMEM/blob/main/SourceSpace/Plots/MMNmean_ConPatDrugInt_Source.R):

	source("MMNmean_ConPatDrugInt_Source.R")
  
	![LFP_MMN3_ConPatDrugInt_wfig_RAUD](https://user-images.githubusercontent.com/23748735/178841262-d03f9e92-bef5-4874-8521-69422e209aa7.png)

<br />

*And the principal finding*

2. [Responses to drug (in auditory cortex) in patients are conditional on GABA concentrations](https://github.com/AlistairPerry/FTLDMEGMEM/blob/main/SourceSpace/Plots/LMM_druggabaint.R):

	source("LMM_druggabaint.R")
  
	![LFP_DrugGABAInt_RAUD_scat](https://user-images.githubusercontent.com/23748735/178842801-dfd39e24-2381-40a8-987c-3c1c678fad3b.png)


<br />

## RMarkdown

Included also is [RMarkdown](https://github.com/AlistairPerry/FTLDMEGMEM/blob/main/SourceSpace/Plots/DrugGABAInteraction.Rmd) file and [resultant pdf](https://github.com/AlistairPerry/FTLDMEGMEM/blob/main/SourceSpace/Plots/DrugGABAInteraction.pdf) output, which detail the workflow steps to reproduce the two main findings above. An example:

<br />

<img width="617" alt="Screenshot 2022-07-14 at 12 39 07" src="https://user-images.githubusercontent.com/23748735/178974380-7fe84e0b-df05-4be5-827c-990f6c7f5d14.png">

## Python

Contained are principal findings reproduced with [Python](https://github.com/AlistairPerry/FTLDMEGMEM/blob/main/SourceSpace/Plots/LMM_druggabaint_py.py) Libraries (i.e. pandas, scypy, seaborn, sklearn)

<br />

The added bonus is that using machine learning tools it assesses how well we can predict drug responses based on patients GABA concentrations.

<br />

With Random Forest Regression and Leave-One-Out-Cross-Validation, indeed we find a relationship between __predicted drug responses__ (y-axis) from our model based on __patients actual responses__ (x-axis):

<br />

![DrugGABAInt_LOOCV_RAUD](https://user-images.githubusercontent.com/23748735/180659291-b0b44552-d813-4743-9715-433847650694.png)

<br />

The __Leave-One-Out Cross-Validation (or LOOCV)__ works by fitting a (_Random Forest_) Regression model to all but one subject, and then repeated with each subject left out. By design it assesses how good our model is based upon unseen data (i.e. left out subject).

<br />

We can run this with the following code:

```python
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestRegressor
    
#Create LOOCV procedure
CrossVal = LeaveOneOut()
    
#Create model
RF_model = RandomForestRegressor(random_state=1)
    
#Evaluate model
CrossVal_scores = cross_val_score(RF_model, X, y, scoring='neg_mean_absolute_error', cv=CrossVal)
```
* _Where X and y represent independent and dependent variables, respectively._

<br />

And we can use the mean absolute error (MAE) to assess model performance which is appropriate for regression:

```python
#Report performance
print('Mean absolute error is: %.3f (%.3f)' % (np.mean(CrossVal_scores), np.std(CrossVal_scores)))
```

<br />


Which reveals the mean error across cross-validations in predicting patients drug responses is __0.033__ (with SD = 0.024).
