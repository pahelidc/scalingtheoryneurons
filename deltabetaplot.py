#DELTA BETA PLOT

#Author: Samhita Mallavarapu,Paheli Desai-Chowdhry
#February-May 2022

#This is a Python Script that produces a plot of the distribution of delta beta (asymmetry level)
#This script takes as its input dat files with a list of beta1 and beta 2 #to be plotted after filtering
#This script takes this data and plots a distribution of the difference scale factor (measure of asymmetry)
#To run this code in the terminal, write python deltabetaplot.py
#The output file is an image with the plot of the distribution

import numpy as np
import scipy as sp
from scipy import stats
import matplotlib.pyplot as plt
import math
import os
import sys

beta2 = np.loadtxt('beta2.dat')
beta1 = np.loadtxt('beta1.dat')
beta = (beta1 - beta2)/2 #calculates the difference scale factor

mean = np.mean(beta)
mean = np.around(mean, decimals =2)#rounds to two decimal places
#mean = np.abs(mean)

std = np.std(beta)
std = np.around(std, decimals =2)#rounds to two decimal places

sem = sp.stats.sem(beta)
sem = np.around(sem, decimals =3)#rounds to three decimal places

plt.plot(1)
plt.hist(beta, bins=20, alpha = 0.75, color = 'salmon', edgecolor = 'k', linewidth = .5) #can change color depending on cell type
plt.title('Motoneurons', fontsize = 27)#change depending on cell type
plt.ylabel('Frequency', fontsize = 23)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xticks(np.arange(-2.5,2.7,0.5)) #modify depending on data 
#plt.xticks(np.arange(min(beta)+0.053, max(beta), 1.0))
plt.xlabel(r'$\Delta \beta$', fontsize = 23)
plt.axvline(mean, color='k', linestyle = 'solid', linewidth=2) #solid line for the mean
plt.text(0.60,590, r'$\mu = %s$' % ('0.00'), fontsize = 25)#move position of text depending on data/visuals
plt.text(0.60,535, r'$\sigma= %s$' % (std), fontsize = 25)
plt.errorbar(mean, 450, xerr= 2*sem, yerr=0, fmt='.k', linewidth=1) #adds error bar for uncertainty measure of the mean

plt.savefig('moto_delta_beta-dec22.png') #change file name by cell type and data
