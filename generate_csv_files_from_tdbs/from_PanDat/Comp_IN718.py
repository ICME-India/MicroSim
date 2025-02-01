# -*- coding: utf-8 -*-
"""
Spyder Editor
@author: Swapnil

"""

import numpy as np
import pandas as pd

#%%
#input dat filename generated from PanDat
fl = open("Comp_IN718.dat","r")

#output filename
f_name = "Comp_IN718.csv"

#The element and phase names should be same as they appear in .dat file from PanDat. 
#The order can be different compared to the .dat file.
#You may choose the order of elements as you want them to appear in the .csv file generated from this code.
#The solvent element is not written out in the .csv file.
elem = ["Si","Co","Ti","Cr","Mn","Fe","Al","Nb","Mo"]
ph = ["Fcc","Liquid"]  #order of output column 

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

num_comp = len(elem)
num_col = num_comp

#%%
#create a Pandas dataframe with empty Hessian and Temp columns
comp = [[[0] * (num_col) for i in range(num_T)] for j in range(num_phases)]
data = [None]*num_phases      
for l in range(num_phases):
    data[l] = pd.DataFrame(comp[l],columns=elem,index=T)
    data[l].index.name = "Temp"


#%%
#populate Hessian data into the dataframe
for count_T in range(num_T):
    for count_phase in range(num_phases):
        count_filled_col = 0
        #while count_filled_col < num_ind_col:
        if count_phase == 0:
            i = pos_T[count_T]+2
        else:
            i = pos_T[count_T]+2 + (num_col+1)*count_phase
        for j in range(num_col):
            #print(line_data[i][0])
            if line_data[i][0] in data[count_phase].columns: # to check if col exists in dataframe
                data[count_phase].loc[T[count_T],line_data[i][0]] = line_data[i][1]
                #print(data[count_phase].loc[T[count_T],line_data[i][0]])
            i += 1       
 


#%%
col_names_out = [None]*(num_col*num_phases)

count = 0
for i in range(num_phases):
    for j in range(num_col):
        #loc = i*(num_comp-1)+j
        col_names_out[count] = str(elem[j])+"@"+str(ph[i])
        #print(col_names_out[count])
        count += 1
               
        
#%%
#combine all phases into one file
data_all = pd.concat([data[phase.index(ph[0])],data[phase.index(ph[1])]], axis=1,join='inner')

    
#%%

data_all.columns = col_names_out
data_all = data_all.sort_values(by="Temp") #sort in ascending order of Temp
data_all.to_csv(f_name)



