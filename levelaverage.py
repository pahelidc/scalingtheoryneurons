#LEVEL AVERAGE

#Author: Paheli Desai-Chowdhry
#May 2023

#This is a Python Script that takes raw data and averages the avg beta and delta beta for each level
#This script takes as its input dat files of the raw data, mloutput
#This script takes this data and loops through, finding the avg asymmetric radius scale factors at each level
#It also outputs a data file of the averaged over level data
#Note that it takes the absolute leaf number as input, since it needs an integer, but then outputs the relative leaf number 
#To run this code in the terminal, write python levelaverage.py  datfile.dat label
#The output file is a data file with the averaged data

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import math
import os
import sys

data = np.loadtxt(sys.argv[1]) #takes raw data as input

avgbeta = data[:,0]
deltabeta = data[:,1]
leafno = data[:,4] #can change depending on where it is in data  
#rellevel = data[:,6]

maxleafno = np.amax(leafno)
print(maxleafno)
n = len(leafno)

avgbetaavg = np.zeros(int(maxleafno)+1)#extra to make sure level 0 is counted
avgbetadiff = np.zeros(int(maxleafno)+1)
level = np.zeros(int(maxleafno)+1)

for j in range(int(maxleafno+1)):
    avgbetaavgtemp = 0
    avgbetadifftemp = 0
    countavg = 0
    for i in range(n):
        if leafno[i] == j:
            avgbetaavgtemp = avgbetaavgtemp + avgbeta[i]
            avgbetadifftemp = avgbetadifftemp + deltabeta[i]
            countavg = countavg + 1#adding to count
            if countavg != 0: 
                avgbetaavg[j] = avgbetaavgtemp/countavg#test
                avgbetadiff[j] = avgbetadifftemp/countavg
                level[j] = np.log2(maxleafno/j)#stores the level for each avg


levelavg  = np.c_[avgbetaavg, avgbetadiff, level]

levelavg = levelavg[~(levelavg==0).all(1)] #removes the zero values from the array

np.savetxt('levelavg_%s.dat'% sys.argv[2], levelavg)#label for file as second input
np.savetxt('test_%s.dat'% sys.argv[2], leafno)#label for file as second input
