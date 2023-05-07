#PIXEL-TO-BRANCHING-RATIOS

#Author: Paheli Desai-Chowdhry
#January 2019-April 2023

#This is a Python Script that extracts scaling ratios of neuron branches as well as asymmetric scale factors
#This script takes as its input a standard neuron reconstruction in the swc format
#This standard reconstruction has a list of pixels with coordinates and radius
#This script takes this data and organizes it in branches with branch and parent labels
#From this, it calculates the radius and length scaling ratios
#To run this code in the terminal, write python swc-to-ratios.py filename.swc label group(classification)
#The classification is to prepare data that can be used for cell type classification machine learning methods 
#The label is so that you can use this to run on multiple swc files in the same directory
#The output files for scaling ratios will have a respective number labelling, such as a number labelling

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import math
import os
import sys


#This block takes the file and counts the number of numerical entries, or pixels

numlines = 0

with open(sys.argv[1]) as f:
    for line in f:
        if line != '/n' and not line.startswith('#') and line.strip():
            numlines += 1

print (numlines) #tells you how long the file is 

data = np.loadtxt(sys.argv[1])

pix_id = data[0:numlines, 0] #extracts the pixel ID number from the file
art_id = data [0:numlines, 1] #extracts the identifier/type where 1 is soma, 2 is axon, 3 is a basal dendrite, and 4 is an apical dendrite
x = data [0:numlines, 2] #extracts the x spatial coordinate
y = data [0:numlines, 3] #extracts the y spatial coordinate
z = data [0:numlines, 4] #extracts the z spatial coordinate
radius = data [0:numlines, 5] #extracts the radius at that coordinate
pix_par_id = data [0:numlines, 6] #extracts the parent pixel

np.savetxt('pixidtest.dat', pix_id)
np.savetxt('pix_par_idtest.dat', pix_par_id)

pix_par_id[0] = 0 #sets the first parent pixel to be 0

#Extracts the pixel IDs where branches occur

branching_id = np.zeros(numlines)
np.savetxt('branchidtest.dat', branching_id)

count1 = 0

for i in range (numlines):
    if (pix_id[i] - pix_par_id[i]) >= 2: #Find where there is a gap between pixel and parent
        branching_id[count1] = branching_id[count1] + pix_id[i]
    if branching_id[count1] != 0:
        count1 = count1 + 1

#Here, we find the child pixel IDs where the branchings occur and the pixel IDs of the parent IDs
#each parent ID is associated with at least two separate child IDs

pix_child_id = np.zeros(numlines)

pix_child_parent_id = np.zeros(numlines)

count2 = 0

for j in range(count1): 
    for i in range(numlines):
        if pix_par_id[int(branching_id[j])-1] == pix_par_id[i]: #Finds the parents of where the branches occur
            pix_child_parent_id[count2] = pix_par_id[i]
            pix_child_id[count2] = i+1 #Finds where the children of the parents occur
            count2 = count2 + 1

#separates the branches 
vessel_id = np.zeros(numlines, dtype=int)

vessel_id[0] = 1
vessel_id[1] = 1

count3 = 1

for i in range(2,numlines):
    for j in range(count2):
        if pix_id[i] == pix_child_id[j]: #If this is where one of the branchinds occur
            count3 = count3 + 1
            break
    vessel_id[i] = count3

vessel_parent_id_list = np.zeros(numlines, dtype=int)

#Extraction of parent ids


for j in range(numlines):
    if vessel_id[j] == 1:
      vessel_parent_id_list[j] = 1
    for i in range(count2):	
        if pix_id[j] == pix_par_id[int(pix_child_id[i])-1]: #Finds where the parents of the branches occur
            vessel_parent_id_list[int(pix_child_id[i])-1] = vessel_id[j] 
            

for k in range(numlines):
    if vessel_parent_id_list[k]==0:
        for j in range(k):
            if pix_par_id[k] == j+1: #For the remaining parent IDs not identified by branches, finds the parent pixels and labels the vessel parent IDs as the same, as they adjacent and in the same vessel
                vessel_parent_id_list[k] = vessel_parent_id_list[j] 
              
               

#Breaks up the parents vessels into levels

vessel_parent_id = np.zeros(count3)

vessel_parent_id[0] = 1
count4 = 1

for j in range(1,numlines):
    if count4 == count3:
       break
    vessel_parent_id[count4] = vessel_parent_id_list[j+1] #labels parent vessels
    if vessel_id[j] != vessel_id[j+1]: #if these is a break in vessels
        count4 = count4+1	    

#Creates the output file, with 4 columns, organized by branch/"vessel"

swc_output = np.zeros((count3, 5))

#Just labels in order
for i in range(count3):
    swc_output[i,0] = i +1
    
#Parents
for i in range(count3):
    swc_output[i,1] = vessel_parent_id[i] #Labels parents in the output so the ratio can be found
    
    
for j in range(count3):
    radius_avg_j = 0
    count4_j = 0
    for i in range(numlines):
        if vessel_id[i] == j+1: 
            count4_j = count4_j +1
            radius_avg_j = radius_avg_j + radius[i] #Sums up randius in a branch
    if count4_j != 0:
        radius_avg_j = radius_avg_j/count4_j #Divides by the total number to find the average
    if count4_j == 0:
        radius_avg_j = radius_avg_j/1
    swc_output[j,2] = swc_output[j,2] + radius_avg_j #Finds the average radius for the output in each branch
    
for j in range(count3):
    length_j = 0
    for i in range(numlines-1):
        if vessel_id[i] == j+1 and vessel_id[i+1] == j+1: 
            length_j = length_j + math.sqrt((x[i+1]-x[i])**2 + (y[i+1]-y[i])**2 + (z[i+1]-z[i])**2) #Finds length by the euclidean distance between pixels
    swc_output[j,3] = swc_output[j,3] + length_j

# distal branches #

num_branches = count3

num_chil = np.zeros(num_branches)
num_distal = np.zeros(num_branches)

for i in range(num_branches):
    for j in range(num_branches):
        if swc_output[i,0] == vessel_parent_id[j] and i != j:#swc_output[i,0] = vessel_id[i]
            num_chil[i] = num_chil[i] + 1


for k in range(num_branches):
    if num_chil[k] == 0:
        num_distal[int(vessel_parent_id[k])-1] = num_distal[int(vessel_parent_id[k])-1] + 1


for m in range(1,num_branches):
    for k in range(num_branches):
        if num_distal[k] == m and vessel_parent_id[k] != swc_output[k,0]:#swc_output[k,0] = vessel_id[k]
            num_distal[int(vessel_parent_id[k])-1] = num_distal[int(vessel_parent_id[k])-1] + num_distal[k]

for i in range(count3):
    swc_output[i,4] = num_distal[i]

#Getting the ratios

ratios = np.zeros((count3, 3))

#Print in first column the parent branch of the ratio, ie, n+1 

for i in range(count3):
    ratios[i,0] = i+1
    for j in range(count3):
        if swc_output[i,0] == swc_output[j,1] and swc_output[i,2] != 0 and swc_output[i,3]!= 0: #j denotes the daughter branch; checking to see if its parent i exists and finding the corresponding daughter to parent ratio 
            ratios[i,1] = swc_output[j,2]/swc_output[i,2] #radius ratio
            ratios[i,2] = swc_output[j,3]/swc_output[i,3] #length ratio
            
count4 = 0
for i in range(count3):
    if ratios[i,1] != 0 and ratios[i,2] != 0: #total number of branches with nonzero ratios
        count4 = count4 + 1

#If at any point the code does not work, can remove comment to save files that have printed intermediate steps 

#np.savetxt('pix_id.dat', pix_id)
#np.savetxt('pix_par_id.dat', pix_par_id)
#np.savetxt('pix_child_id.dat', pix_child_id)
#np.savetxt('pix_child_parent_id.dat', pix_child_parent_id)
#np.savetxt('branching_id.dat', branching_id)
#np.savetxt('vessel_id.dat', vessel_id)
#np.savetxt('vessel_parent_id_list.dat', vessel_parent_id_list)
#np.savetxt('vessel_parent_id.dat', vessel_parent_id)
#np.savetxt('ratios.dat', ratios)

#Intermediate output: 

int_output = np.c_[vessel_id, vessel_parent_id_list, x, y, z, radius]
#np.savetxt('int_output_%s.dat'% sys.argv[2], int_output)

#Final output: Angicart-style output (arranged by branch) and list of ratios 

np.savetxt('swc_output_%s.dat'% sys.argv[2], swc_output)
#np.savetxt('radius_ratio_list_%s.dat'% sys.argv[2], radius_ratio_list)
#np.savetxt('length_ratio_list_%s.dat'% sys.argv[2], length_ratio_list)

#ML cell-type classification data processing 

#We want to extract the asymmetric average and difference scale factors

beta_avg = np.zeros(count3)#average radius scale factor
gamma_avg = np.zeros(count3)#average length scale factor
beta_diff = np.zeros(count3)#difference radius scale factor
gamma_diff = np.zeros(count3)#difference length scale factor
leaf_number_beta = np.zeros(count3) #number of distal tips
leaf_number_gamma = np.zeros(count3) #number of distal tips
asym_ratio = np.zeros(count3) #ratio of leaf numbers of daughters 

beta_1 = np.zeros(count3)#Use these to calculate scaling exponents (after filtering)
beta_2 = np.zeros(count3)
gamma_1 = np.zeros(count3)
gamma_2 = np.zeros(count3)

for i in range(count3):
    for j in range(count3):
        if vessel_parent_id[i] == vessel_parent_id[j] and i != j and swc_output[int(vessel_parent_id[i])-1,2]!= 0 and swc_output[int(vessel_parent_id[i])-1,3]!= 0:
            beta_avg[i] = (swc_output[i,2]+swc_output[j,2])/(2*swc_output[int(vessel_parent_id[i])-1,2])
            beta_1[i] = (swc_output[i,2])/(swc_output[int(vessel_parent_id[i])-1,2])#test for cons based exponents           
            beta_2[i] = (swc_output[j,2])/(swc_output[int(vessel_parent_id[i])-1,2])#test for cons based exponents
            leaf_number_beta[i] = swc_output[int(vessel_parent_id[i])-1,4]
            gamma_avg[i] = (swc_output[i,3]+swc_output[j,3])/(2*swc_output[int(vessel_parent_id[i])-1,3])
            gamma_1[i] = (swc_output[i,3])/(swc_output[int(vessel_parent_id[i])-1,3])#test for cons based exponents
            gamma_2[i] = (swc_output[j,3])/(swc_output[int(vessel_parent_id[i])-1,3])#test for cons based exponents
            leaf_number_gamma[i] = swc_output[int(vessel_parent_id[i])-1,4]
            beta_diff[i] = (swc_output[i,2]-swc_output[j,2])/(2*swc_output[int(vessel_parent_id[i])-1,2])
            gamma_diff[i] = (swc_output[i,3]-swc_output[j,3])/(2*swc_output[int(vessel_parent_id[i])-1,3])
            if gamma_diff[i] < 0: #making all the signs of length consistent; radius changes in sign depending on direction of asymmetry
                beta_diff[i] = (swc_output[j,2]-swc_output[i,2])/(2*swc_output[int(vessel_parent_id[i])-1,2])
                gamma_diff[i] = (swc_output[j,3]-swc_output[i,3])/(2*swc_output[int(vessel_parent_id[i])-1,3])

#Asymmetry ratios

for i in range(count3): 
    for j in range(count3):
        if vessel_parent_id[i] == vessel_parent_id[j] and i != j:
            if swc_output[i,4] == 0:
                swc_output[i,4] = 1
            if swc_output[j,4] == 0:
                swc_output[j,4] = 1
            if swc_output[i,4] <= swc_output[j,4]:
                asym_ratio[i] = swc_output[i,4]/swc_output[j,4]
            if swc_output[j,4] <= swc_output[i,4]:
                asym_ratio[i] = swc_output[j,4]/swc_output[i,4]

#This loop below removes those that are identical, i.e.: daughters of the same parent
for i in range(count3):
    for j in range(count3):
        if i != j and beta_avg[i] == beta_avg[j] and gamma_avg[i] == gamma_avg[j] and beta_diff[i] == beta_diff[j] and gamma_diff[i] == gamma_diff[j]:
            beta_avg[j] = 0
            gamma_avg[j] = 0 
            beta_diff[j] = 0
            gamma_diff[j] = 0
            leaf_number_beta[j] = 0
            asym_ratio[j] = 0

#Removes values greater than or equal to 1 (used 0.999 to remove the ratios that are approximately 1 due to the resolution limit of the images)  
for j in range(count3):
    if beta_avg[j] >= 0.999 or beta_1[j] >= 0.999 or beta_2[j] >= 0.999:
        beta_avg[j] = 0
        gamma_avg[j] = 0              
        beta_diff[j] = 0
        gamma_diff[j] = 0
        leaf_number_beta[j] = 0
        asym_ratio[j] = 0

ml_output = np.c_[beta_avg, beta_diff, gamma_avg, gamma_diff, leaf_number_beta, asym_ratio]

ml_output = ml_output[~(ml_output==0).all(1)] #removes the zero values from the array

n = len(ml_output) #extracts the length of the output

#This is specific to the purk/moto classification, adds either 1 or 0 for group label. Can change depending on comparison/data type
if sys.argv[3] == "purk": 
    G = np.ones(n)
    ml_output = np.c_[ml_output, G]
if sys.argv[3] == "moto": 
    G = np.zeros(n)
    ml_output = np.c_[ml_output, G]


#additional data: cummulative for each cell
maxleafno = np.amax(ml_output[:,4])

summary = np.c_[np.mean(ml_output[:,0]), np.mean(ml_output[:,1]), np.amax(ml_output[:,4]), np.mean(ml_output[:,4])] #summary of cell, including maximum and average leaf number

np.savetxt('cellsummary_%s.dat'% sys.argv[2], summary)

#additional data: relative level based on leaf number, based on the approximation that there are 2^(level) leaves at each level

relativelevel = np.zeros(n)#creating a column of the same length as ml output with the relative level 

for i in range(n):
    relativelevel[i] = np.log2(maxleafno/ml_output[i,4])

ml_output = np.c_[ml_output, relativelevel]
np.savetxt('ml_output_%s.dat'% sys.argv[2], ml_output)#saving ml output file with the relative level 

#Data averaged across each level

avgbetaavg = np.zeros(int(maxleafno)+1)#extra to make sure level 0 is counted
avgbetadiff = np.zeros(int(maxleafno)+1)
level = np.zeros(int(maxleafno)+1)

for j in range(int(maxleafno+1)):
    avgbetaavgtemp = 0
    avgbetadifftemp = 0
    countavg = 0
    for i in range(n):
        if ml_output[i,4] == j:
            avgbetaavgtemp = avgbetaavgtemp + ml_output[i,0]
            avgbetadifftemp = avgbetadifftemp + ml_output[i,1] 
            countavg = countavg + 1#adding to count 
            if countavg != 0:	
                avgbetaavg[j] = avgbetaavgtemp/countavg#test
                avgbetadiff[j] = avgbetadifftemp/countavg
                level[j] = np.log2(maxleafno/j)#stores the level for each avg
	

levelavg  = np.c_[avgbetaavg, avgbetadiff, level]

levelavg = levelavg[~(levelavg==0).all(1)] #removes the zero values from the array

np.savetxt('levelavg_%s.dat'% sys.argv[2], levelavg)

#test for scaling exponents, filtering data

for i in range(len(beta_1)):
    if beta_2[i] == 0:
        beta_1[i] = 0

for i in range(len(beta_2)):
    if beta_1[i] == 0:
        beta_2[i] = 0

leaf_number_beta = leaf_number_beta[beta_1 != 0]
beta_1 = beta_1[beta_1 != 0]
beta_2 = beta_2[beta_2 != 0]


np.savetxt('beta1_%s.dat'% sys.argv[2], beta_1)#saves with the cell-type label 
np.savetxt('beta2_%s.dat'% sys.argv[2], beta_2)
np.savetxt('leaf_number_beta_%s.dat'% sys.argv[2], leaf_number_beta)

for i in range(len(gamma_1)):
    if gamma_2[i] == 0: 
        gamma_1[i] = 0

for i in range(len(gamma_2)):
    if gamma_1[i] == 0: 
        gamma_2[i] = 0

leaf_number_gamma = leaf_number_gamma[gamma_1 != 0]

gamma_1 = gamma_1[gamma_1 != 0]
gamma_2 = gamma_2[gamma_2 != 0]

np.savetxt('gamma1_%s.dat'% sys.argv[2], gamma_1)
np.savetxt('gamma2_%s.dat'% sys.argv[2], gamma_2)
np.savetxt('leaf_number_gamma_%s.dat'% sys.argv[2], leaf_number_gamma)
