# -*- coding: utf-8 -*-
"""
Spyder Editor
@author: Swapnil

"""

import numpy as np
import pandas as pd

#%%
#input dat filename generated from PanDat
fl = open("Hessian_IN718.dat","r")
system_name = "IN718" # for output filename

#The element and phase names should be same as they appear in .dat file from PanDat. 
#The order can be different compared to the .dat file.
#You may choose the order of elements as you want them to appear in the .csv file generated from this code.
#The solvent element is not written out in the .csv file.
elem = ["Si","Co","Ti","Cr","Mn","Fe","Al","Nb","Mo"]



#%%
#create col names based on user-defined elem order
num_comp = len(elem)

num_col = int(num_comp**2)
num_ind_col = int(num_col - (num_comp-1)*(num_comp)/2)

col_names_out = [None]*num_ind_col
count = 0
for i in range(num_comp):
    for j in range(num_comp):
        #loc = i*(num_comp-1)+j
        if i==j:
            col_names_out[count] = "HSN["+str(elem[i])+","+str(elem[j])+"]"
            count += 1
for i in range(num_comp):
    for j in range(num_comp):    
        if i<j:
            col_names_out[count] = "HSN["+str(elem[i])+","+str(elem[j])+"]"
            count += 1
    

#%%
data = []
for i in fl.readlines():
  data.append(i) 
  ## each line is stored as an element in data list


#%%
#find location, value and num of T
# find num of phases
line_data = [None]*len(data)
T     = []
pos_T = []
phase = []
col_names = []

num_T = 0
num_phases = 0

for i in range(len(data)):
    line_data[i] = data[i].split()
    if len(line_data[i])!=0:
        if line_data[i][0]=="Temperature":
            pos_T.append(i)
            T.append(line_data[i][1])
            num_T += 1 
        if num_T == 1 and line_data[i][0] == "Phase":
            phase.append(line_data[i][3])
            num_phases += 1
    i += 1  
        
#%%
#create a Pandas dataframe with empty Hessian and Temp columns
hsn = [[[0] * (num_ind_col) for i in range(num_T)] for j in range(num_phases)]
data = [None]*num_phases      
for l in range(num_phases):
    data[l] = pd.DataFrame(hsn[l],columns=col_names_out,index=T)
    data[l].index.name = "Temp"


#%%
#populate Hessian data into the dataframe
for count_T in range(num_T):
    for count_phase in range(num_phases):
        count_filled_col = 0
        #while count_filled_col < num_ind_col:
        if count_phase == 0:
            i = pos_T[count_T]+3
        else:
            i = pos_T[count_T]+2 + (num_col+7)*count_phase
        for j in range(num_col):
            if line_data[i][0] in data[count_phase].columns: # to check if col exists in dataframe
                data[count_phase].loc[T[count_T],line_data[i][0]] = line_data[i][1]
                #print(data[count_phase].loc[T[count_T],line_data[i][0]])
            i += 1       
  
#%%
#modify column names
col_names_out_mod = [[0] * (num_ind_col) for i in range(num_phases)]

for i in range(num_phases):
    for j in range(num_ind_col):
        col_names_out_mod[i][j] = str(col_names_out[j]) + "@"+phase[i]
    data[i].columns = col_names_out_mod[i]
              
#%%
f_name = [None]*num_phases

for l in range(num_phases):
    f_name[l] = "HSN_"+str(phase[l])+"_"+str(system_name)+".csv"
    data[l] = data[l].sort_values(by="Temp") #sort in ascending order of Temp
    data[l].to_csv(f_name[l])



