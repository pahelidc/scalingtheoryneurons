#LEAF NUMBER PLOT

#Author: Samhita Mallavarapu,Paheli Desai-Chowdhry
#February-May 2022

#This is a Python Script that produces a scatter plot of the relationship between leaf number and asymmetry
#This script takes as its input dat files with a list of beta1 and beta 2 and leaf number, which should all be the same length
#This script takes this data and plots a scatter plot with leaf number on the x axis and the absolute value of delta beta in the y-axis
#These plots show the relationship between magnitude of asymmetry and relative position in the cell
#To run this code in the terminal, write python leafnumberplot.py
#The output file is an image with the plot of the distribution

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

beta1 = np.loadtxt('beta1.dat')#change depending on file names of data
beta2 = np.loadtxt('beta2.dat')

betadiff = np.subtract(beta1,beta2)/2
std = np.std(betadiff)#calculates the standard deviation as reported in the delta beta plots

betadiff = np.abs(betadiff) #here, we are focusing on the magnitude of asymmetry 

leafnobeta = np.loadtxt('leaf_number_beta_moto.dat')#change depending on file names of data

plt.figure(3)
plt.title('Motoneurons', fontsize = 27)#change depending on cell type
plt.ylabel(r'$|\Delta \beta|$', fontsize = 25)
plt.xlabel('Leaf Number', fontsize = 25)
plt.xlim(left = -50, right = 1250)
plt.ylim(bottom = -0.15, top = 3)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.axhline(2*std, color='k', linestyle = 'dashed', linewidth=2, label = 'Asymmetry Line') #plots a dashed line to show the point that is two standard deviations away from symmetry 
plt.scatter(leafnobeta.astype(int), betadiff, color = 'salmon', linewidths = 2,
            marker ="^",
            edgecolor ="maroon",
            s = 100)
plt.text(750, 0.45, 'Asymmetry Line', fontsize = 20) #change position depending on data
plt.savefig('leafnobeta_moto-dec22.png')#change name depending on cell type and data


