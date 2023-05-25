#SCATTER BETA VS LEVEL

#Author: Paheli Desai-Chowdhry
#April-May 2023

#This is a Python Script that produces a scatter plot of the relationships between relative leaf number and the asymmetric scaling ratios
#This script takes as its input dat files the level avg data
#This script takes this data and plots a scatter plot with leaf number on the x axis and the asymmetric scaling ratios on the y axes
#To run this code in the terminal, write python scatterbetavslevel.py datfile1.dat datfile2.dat.....etc label
#Can be adjusted for different data and additional groups; colors and labels on the plots can be changed
#The output files are 2 images, for avg and delta beta

from mpl_toolkits import mplot3d
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import math
import os
import sys
import csv

group1 = np.loadtxt(sys.argv[1])

beta_avg_group1 = group1[:,0]
beta_diff_group1 = np.abs(group1[:,1])
leafno_group1 = group1[:,2]#relativeleafno

group2 = np.loadtxt(sys.argv[2])

beta_avg_group2 = group2[:,0]
beta_diff_group2 = np.abs(group2[:,1])
leafno_group2 = group2[:,2]#relativeleafno

group3 = np.loadtxt(sys.argv[3])

beta_avg_group3 = group3[:,0]
beta_diff_group3 = np.abs(group3[:,1])
leafno_group3 = group3[:,2]#relativeleafno

group4 = np.loadtxt(sys.argv[4])

beta_avg_group4 = group4[:,0]
beta_diff_group4 = np.abs(group4[:,1])
leafno_group4 = group4[:,2]

group5 = np.loadtxt(sys.argv[5])

beta_avg_group5 = group5[:,0]
beta_diff_group5 = np.abs(group5[:,1])
leafno_group5 = group5[:,2]

group6 = np.loadtxt(sys.argv[6])

beta_avg_group6 = group6[:,0]
beta_diff_group6 = np.abs(group6[:,1])
leafno_group6 = group6[:,2]

#group7 = np.loadtxt(sys.argv[7])

#beta_avg_group7 = group7[:,0]
#beta_diff_group7 = np.abs(group7[:,1])
#leafno_group7 = group7[:,2]

#group8 = np.loadtxt(sys.argv[8])

#beta_avg_group8 = group8[:,0]
#beta_diff_group8 = np.abs(group8[:,1])
#leafno_group8 = group8[:,2]

#group9 = np.loadtxt(sys.argv[9])

#beta_avg_group9 = group9[:,0]
#beta_diff_group9 = np.abs(group9[:,1])
#leafno_group9 = group9[:,2]

#group10 = np.loadtxt(sys.argv[10])

#beta_avg_group10 = group10[:,0]
#beta_diff_group10 = group10[:,1]
#leafno_group10 = group10[:,2]

plt.figure(1)
plt.scatter(leafno_group1, beta_avg_group1, color = "red", label="axons", marker ="^")
plt.scatter(leafno_group2, beta_avg_group2, color = "gold", label= "microglia", marker ="^")
plt.scatter(leafno_group3, beta_avg_group3, color = "green", label="motoneurons", marker ="^")
plt.scatter(leafno_group4, beta_avg_group4, color = "deepskyblue", label="purkinje", marker ="^")
plt.scatter(leafno_group5, beta_avg_group5, color = "blue", label = "MSN", marker ="^")
plt.scatter(leafno_group6, beta_avg_group6, color = "magenta", label = "Pyramidal", marker ="^")
#plt.scatter(leafno_group7, beta_avg_group7, color = "blue")
#plt.scatter(leafno_group8, beta_avg_group8, color = "darkviolet", label = "Pyramidal")
#plt.scatter(leafno_group9, beta_avg_group9, color = "darkviolet")
#plt.scatter(leafno_group10, beta_avg_group10, color = "darkviolet")
plt.title('Radius Average Scale Factor vs Level', fontsize = 15)
plt.xlabel('Relative Leaf Number (Level)', fontsize = 15)
plt.ylabel(r'Radial Average Scale Factor ($\bar \beta$)', fontsize = 15)
plt.xlim(right = 12)
plt.legend(loc ="lower right")
plt.savefig('AvgBetavslevel_%s.png'% sys.argv[7])

plt.figure(2)
plt.scatter(leafno_group1, beta_diff_group1, color = "red", label="axons", marker ="^")
plt.scatter(leafno_group2, beta_diff_group2, color = "gold", label= "microglia", marker ="^")
plt.scatter(leafno_group3, beta_diff_group3, color = "green", label="motoneurons", marker ="^")
plt.scatter(leafno_group4, beta_diff_group4, color = "deepskyblue", label="purkinje", marker ="^")
plt.scatter(leafno_group5, beta_diff_group5, color = "blue", label = "MSN", marker ="^")
plt.scatter(leafno_group6, beta_diff_group6, color = "darkviolet", label = "Pyramidal", marker ="^")
#plt.scatter(leafno_group7, beta_avg_group7, color = "blue")
#plt.scatter(leafno_group8, beta_avg_group8, color = "darkviolet", label = "Pyramidal")
#plt.scatter(leafno_group9, beta_avg_group9, color = "darkviolet")
#plt.scatter(leafno_group10, beta_avg_group10, color = "darkviolet")
plt.title('Radius Difference Scale Factor vs Level', fontsize = 15)
plt.xlabel('Relative Leaf Number (Level)', fontsize = 15)
plt.ylabel(r'Radial Difference Scale Factor ($|\Delta \beta|$)', fontsize = 15)
plt.xlim(right = 12)
plt.legend(loc ="upper right")
plt.savefig('DeltaBetavslevel_%s.png'% sys.argv[7])
