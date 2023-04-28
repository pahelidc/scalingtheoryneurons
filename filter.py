#FILTERING

#Author: Samhita Mallavarapu,Paheli Desai-Chowdhry
#February 2022

#This is a Python Script that filters the beta1, beta2, and leaf number data to prepare it for analaysis
#This script takes as its input dat files with a list of beta1/beta2 and leaf number, which must all correspond and be of the same length
#This script takes this data and filters it to remove data that is not usable due to the resolution limit of the image
#To run this code in the terminal, write python filter.py
#The output files are the filtered data, which are now ready to process further

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import math
import os
import sys
import csv
from scipy.optimize import fsolve
from scipy import stats

beta1 = np.loadtxt('beta1.dat')#replace these with the specific file names of each cell type 
beta2 = np.loadtxt('beta2.dat')
leafno = np.loadtxt('leaf_number_beta.dat')

n = len(beta1)

for i in range(n):
    if beta1[i] == 1 or beta2[i] == 1:
        beta1[i] = 1
        beta2[i] = 1
        leafno[i] = 0 

beta1 = beta1[beta1 != 1] # removes all these points
beta2 = beta2[beta2 != 1]
leafno = leafno[leafno != 0]


np.savetxt('beta1.dat', beta1moto)
np.savetxt('beta2.dat', beta2moto)
np.savetxt('leaf_number_beta.dat', leafnomoto)



