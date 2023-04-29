#CONFIDENCE INTERVALS

#Author: Paheli Desai-Chowdhry
#March 2023

#This is a Python Script that calculates the confidence intervals for the median scaling exponent  
#This script takes as its input dat files for the scaling exponents as calculated from beta1 and beta2 numerically
#This script takes this data and finds an uncertainty measure for the median of the data
#To run this code in the terminal, write python confinterval.py
#The output shows the high and low limits of the range of the interval, printed 

import numpy as np
import scipy as sp
from scipy import stats
import matplotlib.pyplot as plt
import math
import os
import sys
import csv
from scipy.optimize import fsolve
from scipy.stats import bootstrap

pmoto = np.loadtxt('pmoto.dat')#this is the scaling exponent data as calculated numerically; change based on file name

pmoto = pmoto[pmoto <2]#filter based on the way the data is filtered to remove outliers
pmoto = pmoto[pmoto > 0]

median = np.median(pmoto)

print(median)

pmoto = (pmoto,) #orders the list of data

confint = bootstrap(pmoto, np.median, confidence_level=0.95) #95% confidence interval (can change)


print (confint.confidence_interval) 
