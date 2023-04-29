#SPLIT ASYMMETRIC SCALING EXPONENTS

#Author: Paheli Desai-Chowdhry
#April 2021-April 2023

#This is a Python Script that splits data into symmetric and asymmetric branching junctions, calculates scaling exponents, and plots their distributions
#This script takes as its input dat files with beta1 and beta2, which must be the same length
#This script takes this data and splits it into symmetric and asymmetric data based on an assigned cutoff
#Using the fsolve function in the scipy library, it numerically solves for the scaling exponents
#It then plots the distribution, including a comparison to theoretical predictions and the median with a confidence interval (which can be modified) 
#The output files are images with the plots of the scaling exponent distributions

import numpy as np
import scipy as sp
from scipy import stats
import matplotlib.pyplot as plt
import math
import os
import sys
import csv
from scipy.optimize import fsolve
from scipy import stats

beta1 = np.loadtxt('beta1.dat')
beta2 = np.loadtxt('beta2.dat')

betadiff = np.subtract(beta1,beta2)

betadiff = np.abs(betadiff) #calculates the difference scale factor magnitufe

#splits the data into symmetric and asymmetric data

beta1sym = beta1[betadiff <0.75] #Can change this value depending on the split of the data into symmetric and asymmetric regimes
beta2sym = beta2[betadiff <0.75]

beta1asym = beta1[betadiff>=0.75]
beta2asym = beta2[betadiff>=0.75]

#saves the split data
np.savetxt('beta1sym.dat', beta1sym) 
np.savetxt('beta2sym.dat', beta2sym)
np.savetxt('beta1asym.dat', beta1asym)
np.savetxt('beta2asym.dat', beta2asym)

#calculates the length of the symmetric and asymmetric data
n = len(beta1sym)
m = len(beta1asym)

psym = np.zeros(n)
pasym = np.zeros(m)

#Numerically solves for the scaling exponent

for i in range(n):
  func_i = lambda psym_i: 1 - beta1sym[i]**psym_i - beta2sym[i]**psym_i
  psym_i = fsolve(func_i, 0.5)
  psym[i] = psym_i

for i in range(m):
  func_i = lambda pasym_i: 1 - beta1asym[i]**pasym_i - beta2asym[i]**pasym_i
  pasym_i = fsolve(func_i, 0.5)
  pasym[i] = pasym_i

np.savetxt('psym.dat', psym)
np.savetxt('pasym.dat', pasym)

#Filters the data
psym = psym[psym > 0]
pasym = pasym[pasym > 0]

psym = psym[psym < 10]
pasym = pasym[pasym < 2]

#Calculates statistical measures of each dataset
mean = np.mean(psym)
mean = np.around(mean, decimals =3)
median = np.median(psym)
median = np.around(median, decimals =2)
std = np.std(psym)
std = np.around(std, decimals =3)
sem = sp.stats.sem(psym)
sem = np.around(sem, decimals =5)


mean2 = np.mean(pasym)
mean2 = np.around(mean2, decimals =3)
median2 = np.median(pasym)
median2 = np.around(median2, decimals =2)
std2 = np.std(pasym)
std2 = np.around(std2, decimals =3)
sem2 = sp.stats.sem(pasym)
sem2 = np.around(sem2, decimals =5)

ci = (1.43-1.37)/2 #replace this with the range of the confidence interval as calculated byy the confidence interval script
center = (1.43+1.37)/2

plt.rc('xtick',labelsize=13)
plt.rc('ytick',labelsize=13)

#Plots the symmetric data
plt.figure(1)
plt.hist(psym, bins=30 , alpha = 0.25, color = 'blueviolet')
plt.title('Motoneuron Symmetric Radius Scaling Exponents', fontsize = 17.5)
plt.ylabel('Frequency', fontsize = 18)
plt.xlabel('Power', fontsize = 18)
#plt.xlim(right = 5)
#plt.ylim(top = 150)
plt.axvline(median, color='k', linestyle = 'solid', linewidth=2, label = 'Median')
plt.axvline(0.75, color='r', linestyle = 'dashed', linewidth=3, label = r'$\bf{P_{TOT}^*}$')#Shows the predictions based on theory for comparison
plt.axvline(1.50, color='orange', linestyle = 'dashed', linewidth=3, label = r'$\bf{P_{MAX/MIN}^*}$')#Shows the predictions based on theory for comparison
#plt.axvline(2, color='forestgreen', linestyle = 'dashed', linewidth=3, label = r'$\bf{P}$')#Shows the predictions based on theory for comparison
#plt.axvline(2.5, color='b', linestyle = 'dashed', linewidth=3, label = r'$\bf{T, \epsilon = 0}$')#Shows the predictions based on theory for comparison
#plt.axvline(3, color='m', linestyle = 'dashed', linewidth=3, label = r'$\bf{T, \epsilon = \frac{1}{2}}$')#Shows the predictions based on theory for comparison
plt.text(5.9,170, r'$Med = %s$' % (median), fontsize = 25) #chang position depending on data
plt.errorbar(center, 225, xerr= ci, yerr=0, fmt='.k', linewidth=3.5)#change position depending on data
plt.legend(loc= 'upper right', prop={'size': 25})
plt.savefig('moto_sym_cons_rad_exp.png')#change name of output file depending on data/cell type

ci2 = (1.07-0.82)/2#replace this with the range of the confidence interval as calculated byy the confidence interval script
center2 = (1.07+0.82)/2

#Plots the asymmetric data
plt.figure(2)
plt.hist(pasym, bins='auto' , alpha = 0.25, color = 'yellow')
plt.title('Motoneuron Asymmetric Radius Scaling Exponents', fontsize = 17.5)
plt.ylabel('Frequency', fontsize = 18)
plt.xlabel('Power', fontsize = 18)
plt.xlim(right = 4)
plt.ylim(top = 30)
plt.axvline(median2, color='k', linestyle = 'solid', linewidth=2, label = 'Median')
plt.axvline(0.75, color='r', linestyle = 'dashed', linewidth=3, label = r'$\bf{P_{TOT}^*}$')
plt.axvline(1.50, color='orange', linestyle = 'dashed', linewidth=3, label = r'$\bf{P_{MAX/MIN}^*}$')
#plt.axvline(2, color='forestgreen', linestyle = 'dashed', linewidth=3, label = r'$\bf{P}$')
#plt.axvline(2.5, color='b', linestyle = 'dashed', linewidth=3, label = r'$\bf{T, \epsilon = 0}$')
#plt.axvline(3, color='m', linestyle = 'dashed', linewidth=3, label = r'$\bf{T, \epsilon = \frac{1}{2}}$')
#plt.text(6.6, 150, r'$\mu = %s$' % (mean), fontsize = 25)#change position depending on data
plt.text(2.4, 14, r'$Med = 0.90$', fontsize = 25)#change position depending on data
plt.errorbar(center2, 9, xerr= ci2, yerr=0, fmt='.k', linewidth=3.5)
plt.legend(loc= 'upper right', prop={'size': 25})
plt.savefig('moto_asym_cons_rad_exp.png')#change name of output file depending on data/cell type

#Note: For plots, we can uncomment the lines that add show theoretical predictions for additional comparisons 
