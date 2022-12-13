#!/usr/bin/env python3

"""
This python module contains all statistical data analysis utilities to more easily analysing
PlatoSim data products. 
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
from colorama import Back, Fore, Style
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter
import statsmodels.api as sm
from statsmodels.graphics.gofplots import qqplot

#==============================================================#
#                           FUNCTIONS                          #
#==============================================================#


def colortheme(theme):
    # Select theme and PI color
    if theme == 'r': color = ['tomato', 'red', 'orange']
    if theme == 'g': color = ['limegreen', 'forestgreen', 'yellowgreen']
    if theme == 'b': color = ['royalblue', 'darkcyan', 'dodgerblue']
    if theme == 'm': color = ['deeppink', 'darkviolet', 'm']
    if theme == 'y': color = ['khaki', 'gold', 'orange']
    return color



def plot_modelfit(data, lsFit, model, lsModel='OLS', CI=[0.05], alpha=0.1, theme='b', 
                  x='x', y='y', xlab='x', ylab='y', yerr=False):

    # Select theme and PI color
    color = colortheme(theme)
        
    # Set parameter strings 
    reg, pre = 'x', 'y'
    
    # Select model string using R terminology (Linear model is default)
    string = (r'Adj. R$^{2} \approx$ ' + str(round(lsFit.rsquared_adj,2)) + '\n' +
              r'AIC     $\approx$ '    + str(round(lsFit.aic,1)) + '\n' +
              r'BIC     $\approx$ '    + str(round(lsFit.bic,1)) + '\n' +
              r'Cond.  $\approx$ '     + str(round(lsFit.condition_number, 1)) + '\n')
    
    # Fetch fit coefficients
    theta_val, theta_err = [], []
    for i in range(len(lsFit.params)):
        theta_val, theta_err = round(lsFit.params[i],2), round(lsFit.bse[i],2)
        string += fr'$\theta_{i} \approx {theta_val} \pm {theta_err}$' + '\n'
    
    # Select title and model
    title  = r'Fit model: $y = \epsilon + \theta_0$'
    if model == 'y ~ x': title += fr' + $\theta_1 {x}$'   
    if model == 'y ~ x + z': title += fr' + $\theta_1$ {x} + \theta_2 z$'
    if model == 'y ~ x + I(x**2)': title += fr' + $\theta_1 {x} + \theta_2 {x}^2$'
    if model == 'y ~ x + I(x**2) + I(x**3)': title += fr' + $\theta_1 {x} + \theta_2 {x}^2 + \theta_3 {x}^3$'
    if model == 'y ~ x + I(np.sin(x))': title += fr' + $\theta_1 {x} + \theta_2 \sin({x})$'            

    # Predict response variable
    n = 100
    xpredict = np.linspace(data[reg].min(), data[reg].max(), n)
    xpredict = np.linspace(data[reg].min(), data[reg].max(), n)
    
    # Select predictions
    ypredict    = lsFit.predict(exog=dict(x=xpredict))
    predictions = lsFit.get_prediction()
    
    # PLOTTING
    fig, ax = plt.subplots(figsize=(13,8))
    
    # plot data
    if yerr:
        ax.errorbar(data[reg], data[pre], yerr=data[yerr], marker='.', color='grey', ls='', alpha=0.1, zorder=0)
        ax.plot(data[reg], data[pre], 'ko', ms=2, alpha=alpha, label="Data", zorder=1)
    else:
        ax.plot(data[reg], data[pre], 'ko', ms=2, alpha=alpha, label="Data", zorder=1)
    
    # Allow more CI intervals
    for i in range(len(CI)):
        df_predictions = predictions.summary_frame(alpha=CI[i])
        df_predictions.index = data.x.values

        # Plot CI
        if len(data) < 1e3:
            ax.fill_between(df_predictions.index, 
                            df_predictions.mean_ci_lower, 
                            df_predictions.mean_ci_upper, 
                            alpha=0.2, color=CI_color[i+1], zorder=2)
        ax.plot(data[reg], df_predictions.mean_ci_lower, '-', c=color[i+1], lw=2, zorder=2, label=str((1-CI[i])*100)+'% CI')
        ax.plot(data[reg], df_predictions.mean_ci_upper, '-', c=color[i+1], lw=2, zorder=2)
    
    # Plot PI -> only relevant for OLS
    if len(data) < 1e3:
        ax.fill_between(df_predictions.index, 
                        df_predictions.obs_ci_lower, 
                        df_predictions.obs_ci_upper, 
                        alpha=0.2, color=PI_color, zorder=3)
    ax.plot(data[reg], df_predictions.obs_ci_lower, '-', c=color[-1], lw=2, zorder=2, label=str((1-CI[0])*100)+'% PI')
    ax.plot(data[reg], df_predictions.obs_ci_upper, '-', c=color[-1], lw=2, zorder=2)
    
    # Plot best fit and data
    ax.plot(df_predictions['mean'], color=color[0], label='Fit', zorder=3)
        
    # Settings
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title,fontsize=15)
    ax.legend(fontsize=13, title='',fancybox=True, framealpha=0.8, loc='upper right')
    ax.set_xlim(data[reg].iloc[0], data[reg].iloc[-1])
    
    # Add fit box
    #RMS = round(np.sqrt(np.mean( (data[pre]-df_predictions['mean'])**2 )))
    #string += f'RMS = {RMS} ppm'
    #props = dict(boxstyle='round', facecolor=color[0], alpha=0.3)
    #ax.text(1.02, 0.98, string, transform=ax.transAxes, fontsize=14, 
    #        verticalalignment='top', bbox=props)            
    plt.show()
    
    
    
    
    
def plot_residuals(data, lsFit, theme='b', reg='x', alpha=0.1, lsModel='OLS'):
    
    # Choose correct residuals
    if lsModel == 'OLS':  resid = lsFit.resid
    else: resid = lsFit.resid_pearson
    
    color = colortheme(theme)
    fig, ax = plt.subplots(1,2, figsize=(16,5))

    # Plot residuals squared vs. observations
    ax[0].plot(data[reg], resid**2, 'ko', ms=2, alpha=alpha, zorder=1)
    ax[0].grid(True,   color='grey',   ls='-',  lw=0.5, zorder=2)
    ax[0].axhline(y=0, color=color[0], ls='--', lw=2.0, zorder=3)
    ax[0].set_xlabel(reg)
    ax[0].set_ylabel(r"Residuals squared, $\varepsilon^2$")
    ax[0].set_xlim(data[reg].iloc[0], data[reg].iloc[-1])
    
    # Plot residuals vs. plot
    ax[1].plot(lsFit.fittedvalues, resid, 'ko', ms=2, alpha=alpha, zorder=1)
    ax[1].grid(True,   color='grey',    ls='-',  lw=0.5, zorder=2)
    ax[1].axhline(y=0, color=color[0],  ls='--', lw=2.0, zorder=3)
    ax[1].set_xlabel("Predicted reponse")
    ax[1].set_ylabel(r"Residuals, $\epsilon$")
    ax[1].set_xlim(np.min(lsFit.fittedvalues), np.max(lsFit.fittedvalues))
    plt.show()
    
    
    
    
def plot_standardized_residuals(data, lsFit, K, reg='x', lsModel='OLS'):
    
    # Choose correct residuals
    if lsModel == 'OLS':  resid = lsFit.resid
    else: resid = lsFit.resid_pearson
    
    N = len(data[reg])
    s2 = np.sum(resid**2) / (N-K)
    standardizedResiduals = resid / np.sqrt(s2)
    
    # Plot standardized residuals (QQ-plot)
    fig, ax = plt.subplots(figsize=(10,6))
    qqplot(standardizedResiduals, line='45', ax=ax)
    plt.show()
