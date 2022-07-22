# FTLD, MEG, and MEMANTINE

This repository contains processing, analysis, and visualisation scripts that were performed in the following paper:

* Perry, A., Hughes, L., Adams, N., Naessens, M., Murley, A., Rouse, M., ... & Rowe, J. (2022). The neurophysiological effect of NMDA-R antagonism of frontotemporal lobar degeneration is conditional on individual GABA concentration. Available as preprint at https://doi.org/10.21203/rs.3.rs-1609477/v1

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

Contained are principal findings reproduced with Python Libraries (pandas, scypy, seaborn, sklearn)

Plus bonus of performing Machine learning techniques

