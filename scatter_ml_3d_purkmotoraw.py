#SCATTER 3D PLOTS

#Author: Paheli Desai-Chowdhry
#April-May 2023

#This is a Python Script that produces a scatter plot of the feature spaces of the asymmetric radius scaling ratios and relative leaf number
#This script takes as its input dat files of the raw data, mloutput
#This script takes this data and plots a scatter plot with leaf number on the z axis and the asymmetric scaling ratios on the x/y axes
#It also outputs a 2D plot of just the scaling ratios
#To run this code in the terminal, write python scatter_ml_3d_purkmotoraw.py datfile1.dat datfile2.dat label
#Can be adjusted for different data and additional groups; colors and labels on the plots can be changed
#The output files are 2 images, 3D and 2D feature spaces

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
beta_diff_group1 = group1[:,1]
leafno_group1 = group1[:,6]#relativeleafno, raw

group2 = np.loadtxt(sys.argv[2])

beta_avg_group2 = group2[:,0]
beta_diff_group2 = group2[:,1]
leafno_group2 = group2[:,6]#rel leaf number, raw 

#group3 = np.loadtxt(sys.argv[3])

#beta_avg_group3 = group3[:,0]
#beta_diff_group3 = group3[:,1]
#leafno_group3 = group3[:,2]#relativeleafno

#group4 = np.loadtxt(sys.argv[4])

#beta_avg_group4 = group4[:,0]
#beta_diff_group4 = group4[:,1]
#leafno_group4 = group4[:,2]

fig = plt.figure(figsize = (10, 7))
ax = plt.axes(projection ="3d")

plt.figure(1)
ax.scatter3D(beta_diff_group1, beta_avg_group1, leafno_group1, color = "green", label="Motoneurons")
ax.scatter3D(beta_diff_group2, beta_avg_group2, leafno_group2, color = "deepskyblue", label="Purkinje Cells")
#ax.scatter3D(beta_diff_group3, beta_avg_group3, leafno_group3, color = "green", label="Motoneurons")
#ax.scatter3D(beta_diff_group4, beta_avg_group4, leafno_group4, color = "deepskyblue", label="Purkinje Cells")
plt.title('Motoneurons vs Purkinje Cells Radius Asymmetry Factors with Leaf Number', fontsize = 15)
plt.xlabel('Radial Difference Scale Factor', fontsize = 13)
plt.ylabel('Radial Average Scale Factor', fontsize = 13)
ax.set_zlabel('Relative Leaf Number', fontsize = 13)
ax.legend()
plt.savefig('mlplot_3D_leafno_%s.png'% sys.argv[3])


plt.figure(2)
plt.scatter(beta_diff_group1, beta_avg_group1, color = "green", label="Motoneurons")
plt.scatter(beta_diff_group2, beta_avg_group2, color = "deepskyblue", label="Purkinje Cells")
plt.title('Motoneurons vs Purkinje Cells Radius Asymmetry Factors', fontsize = 15)
plt.xlabel('Radial Difference Scale Factor', fontsize = 13)
plt.ylabel('Radial Average Scale Factor', fontsize = 13)
plt.legend(loc ="lower right")
plt.savefig('mlplot_2D_%s.png'% sys.argv[3])
