#BIND

#Author: Paheli Desai-Chowdhry
#February 2020

#This is a Python Script that combines data from a series of dat files
#This script takes as its input dat files that are vertical lists
#This script takes this data and combines it in vertical stacks
#To run this code in the terminal, write python bind.py
#The output file is the combined data file

import pandas as pd
import numpy as np

rad1 = np.loadtxt('radius_ratio_list1.dat') #can change this depending on the file names, number of files to be combined, etc
rad2 = np.loadtxt('radius_ratio_list2.dat')
rad3 = np.loadtxt('radius_ratio_list3.dat')
rad4 = np.loadtxt('radius_ratio_list4.dat')
rad5 = np.loadtxt('radius_ratio_list5.dat')

vert_stack = np.concatenate((rad1, rad2, rad3, rad4, rad5), axis=None)

np.savetxt('combineddata.dat', vert_stack)
