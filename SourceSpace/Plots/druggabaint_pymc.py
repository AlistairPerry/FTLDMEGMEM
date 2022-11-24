#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 12:24:41 2022

@author: alistairperry
"""


'''

Define function

'''


import theano
from theano.printing import Print
import pymc as pm
import numpy as np
import theano.tensor as T
import scipy


def covariance(sigma, rho):
    C = T.fill_diagonal(T.alloc(rho, 2, 2), 1.)
    S = T.diag(sigma)
    M = S.dot(C).dot(S)
    return M


def analyze_standard(data):
    with pm.Model() as model:
        # priors
        sigma = pm.Uniform('sigma', lower=0, upper=1000, shape=2,
                           testval=[0.133434, 5.9304],  # init with mad
                           transform=None)
        rho = pm.Uniform('r', lower=-1, upper=1,
                         testval=-0.2144021,  # init with Spearman's correlation
                         transform=None)

        # print values for debugging
        rho_p = Print('rho')(rho)
        sigma_p = Print('sigma')(sigma)

        cov = pm.Deterministic('cov', covariance(sigma_p, rho_p))
        cov_p = Print('cov')(cov)

        mult_norm = pm.MvNormal('mult_norm', mu=[10.1, 79.],  # set mu to median
                                cov=cov_p, observed=data.T)

    return model


def analyze_robust(data):
    with pm.Model() as model:
        # priors
        mu = pm.Normal('mu', mu=0., tau=0.000001, shape=2,
                       testval=np.array([10.1, 79.]))  # set mu to median
        sigma = pm.Uniform('sigma', lower=0, upper=1000, shape=2,
                           testval=np.array([0.133434, 5.9304]),  # init with mad
                           transform='interval')
        rho = pm.Uniform('r', lower=-1, upper=1,
                         testval=-0.2144021,  # init with Spearman's correlation
                         transform=None)

        # print values for debugging
        rho_p = Print('rho')(rho)
        sigma_p = Print('sigma')(sigma)

        cov = pm.Deterministic('cov', covariance(sigma_p, rho_p))
        num = pm.Exponential('nu_minus_one', lam=1. / 29., testval=1)
        nu = pm.Deterministic('nu', num + 1)
        cov_p = Print('cov')(cov)
        nu_p = Print('nu')(nu)

        mult_norm = pm.MvStudentT('mult_norm', nu=nu_p, mu=mu,
                                  Sigma=cov_p, observed=data.T)

    return model


'''

Run

'''

import os 
import pandas as pd
import numpy as np


#Load and define variables again

#Set data dir

work_dir = "/Users/alistairperry/Documents/Cambridge/Project/Analysis/LFPs/newmaxfilter_icafixes/LFPs_COH_wrad_allMEGch_wfids/C_Plots/DrugAnalysis/MRS"

os.chdir(work_dir)


#Regions to load from text file

base_fname = 'LFPs_MRSIFGgabacorassoc_meanMMN3_DrugDiff_Pats_fullsample_stats_LMMtableforR_RAUD.txt'

MMNtab = pd.read_table(base_fname)


#Define variables for dataset

data = np.array([MMNtab['GABA'], MMNtab['Diff']])


#Run

robust_model = analyze_standard(data)
with robust_model:
    step = pm.Metropolis()
    robust_trace = pm.sample(500, tune=250, step=step, random_seed=21412, progressbar=False)

print(pm.summary(robust_trace))
