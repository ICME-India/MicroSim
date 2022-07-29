
#!/usr/bin/env python
# coding: utf-8

# In[43]:


import sys

import numpy as np
import pandas as pnd
from pandas import DataFrame as df

from pycalphad import Model, Database, calculate, equilibrium , binplot
import pycalphad.variables as v

import sympy as sym
import matplotlib.pyplot as plt

import time
from tinydb import where

sym.init_printing()

from sympy import *

from sympy.codegen.ast import Assignment


# In[68]:


db_sys = Database(sys.argv[1])


# In[69]:


my_phases_sys = list(db_sys.phases.keys())
my_phases_sys = sorted(my_phases_sys) # my_phases_sys contain phases present for the system
print( 'Phases present: ',  my_phases_sys)


# In[70]:


Solid_Phase = sys.argv[2]
Liquid_Phase = sys.argv[3]


# In[71]:


for phas in my_phases_sys:  
 
    if phas == Liquid_Phase:
        my_phases_sys.remove(phas)      
        my_phases_sys.append(phas)      


# In[74]:


my_component_sys = list(db_sys.elements)
for component in my_component_sys:

    if component == '/-':
        my_component_sys.remove(component)              # my_component_sys now doesn't have '/-' as element    


# In[75]:


calc_component_sys = [x for x in my_component_sys if x != 'VA']
calc_component_sys.sort()
print('Components are :',calc_component_sys)


# In[112]:


Compo_calc_1 = sys.argv[4] 
Index_Compo_calc_1 = calc_component_sys.index(Compo_calc_1)
print("Component for Calculation : ",Compo_calc_1)

ind_ex_1 = calc_component_sys.index(Compo_calc_1)    
    

# In[ ]:


Compo_calc_2 = sys.argv[5] 
Index_Compo_calc_2 = calc_component_sys.index(Compo_calc_2)
print("Other Component is : ",Compo_calc_2)

ind_ex_2 = calc_component_sys.index(Compo_calc_2)  

# In[114]:


phases_equi = [Solid_Phase,Liquid_Phase]   # one has to be solid Phase and other has to be Liquid Phase 
phases_calc = sorted(phases_equi)


# In[115]:


for phas in phases_equi:  # phases_equi has Solid phase as 1st item and Liquid as last item 

    if phas == Liquid_Phase:
        phases_equi.remove(phas)
        phases_equi.append(phas)
print("Phases to calculate equilibrium are :", phases_equi[0],phases_equi[1])        

# In[117]:


c_alloy = float(sys.argv[6])

print('Alloy Composition :',c_alloy)





# In[116]:


strt_A = 'd2G_' + Solid_Phase + '_' + Compo_calc_1 + '_' + str(c_alloy)
strt_B = 'd2G_' + Liquid_Phase + '_' + Compo_calc_1 + '_' + str(c_alloy)
strt_C = 'X_' + Compo_calc_1 + '_' + phases_equi[0] + '_' 
strt_D = 'X_' + Compo_calc_1 + '_' + phases_equi[1] + '_' 

#print(strt_A,strt_B,strt_C,strt_D)

Temp_Tol = 0.001
Comp_Tol = 0.0005


# In[118]:


def equi_composition_phase(Temperature,Composition,phase):  
    
    #Function to compute equilibrium phase composition at a Temperature 
    # Inputs to equilibrium Function : Temperature, Composition
    
    eq_comp = equilibrium(db_sys,my_component_sys,my_phases_sys,{v.X(Compo_calc_1):Composition,v.T:Temperature,v.P: 101325},output='GM')
    Phase_at_equi = np.array(eq_comp.Phase.squeeze())
    Phase_at_equi = [x for x in Phase_at_equi if x] # Eliminate Nan Entries
    Phase_at_equi = sorted(Phase_at_equi)
    
    if phase[0] in Phase_at_equi:
    
        phase_compo = np.array(eq_comp.X.where(eq_comp.Phase == phase).sel(component = Compo_calc_1).squeeze())
        phase_compo = phase_compo[~(np.isnan(phase_compo))]
        return phase_compo[0] 
    
    else:
        
        print(phase[0],' : Not an Equilibrium Phase for Composition at the Temperature')
        return


# In[119]:


def phase_in_equilibrium(Temperature,Composition):
    
    # Function that returns equilibrium phases at a Temperature for a given composition 
    
    eq_com = equilibrium(db_sys,my_component_sys,my_phases_sys,{v.X(Compo_calc_1):Composition,v.T:Temperature, v.P: 101325},output='GM')
    
    Phase_at_equi = np.array(eq_com.Phase.squeeze())
    Phase_at_equi = [x for x in Phase_at_equi if x]
    Phase_at_equi = sorted(Phase_at_equi)
        
    return Phase_at_equi


# In[120]:


def Binary_Search(T_Low,T_High,c_alloy): # Binary Search to find One Temp where equilibrium exists 
    
    # T_High : Temperature above which Liq exists 
    # T_Low : Temperature above which Solid exists
    # Between T_High and T_Low , (Liq + Solid) exists found using Binary Search 
    
    # Search stop when Temp_dif acheives tolerance or Equilibria Temp is found(flag = 1)  
    
    flag = 0
    Temp_dif = abs(T_High - T_Low)
    
    T_Current = (T_High + T_Low)*0.5
    phase_equi_cur = phase_in_equilibrium(T_Current,c_alloy)      

    if phase_equi_cur == phases_calc:
        
        flag = 1
    
    while(flag == 0 and Temp_dif >= Temp_Tol):  
        
        phase_equi_cur = phase_in_equilibrium(T_Current,c_alloy)

        if phase_equi_cur == phases_calc:
     
            flag = 1
            break
    
        elif phase_equi_cur == [Liquid_Phase]:
                
            T_High = T_Current
        
        else:
        
             T_Low = T_Current
        
        T_Current = (T_High + T_Low)*0.5
        Temp_dif = abs(T_High - T_Low)     
    
    if(flag == 0 and Temp_dif < Temp_Tol):
        
        print('Equilibrium not Found')
        
        
    return [T_Low,T_Current,T_High,flag]                


# In[121]:


def Binary_Melting(T_Low,T_High,c_alloy,phase_tem):
    
    # Functions to find Liquidus Temp (phase_tem is Liquid)
    # Functions to find Solidus Temp (phase_tem is Solid)
    
    # T_High : phase Tem : Solid,  Temperature at which (Liq + solid) exists 
    # T_Low : phase Tem : Solid,  Temperature at which solid exists 
    
    # T_High : phase Tem : Liquid,  Temperature at which Liquid exists 
    # T_Low : phase Tem : Solid,  Temperature at which (Liq + solid) exists 
    
    # Between T_High and T_Low , (Liq + Solid) exists found using Binary Search till Temp_tol is acheived 
    # Search stop when Temp_dif acheives a tolerance and Equilibria exists at the Temp  found(flag = 1)  
     
    flag = 0
    
    if phase_tem == [Solid_Phase]:
       
        T_Current = T_High
    
    else:
          
        T_Current = T_Low  
    
    phase_equi_cur = phase_in_equilibrium(T_Current,c_alloy)
    
    if phase_equi_cur == phases_calc:
        
        flag = 1
    
    Temp_dif = abs(T_High - T_Low)
    
    while(not(bool(abs(Temp_dif) < Temp_Tol and flag == 1))):  
        
        T_Current = (T_High + T_Low)*0.5
        phase_equi_cur = phase_in_equilibrium(T_Current,c_alloy)
        
        if phase_equi_cur == phases_calc:
            
            flag = 1
            
            if phase_tem == [Solid_Phase]:
                
                T_High = T_Current
                
            else:
                
                T_Low = T_Current
            
        else:
           
            flag = 0
            
            if phase_tem == [Solid_Phase]:
                
                T_Low = T_Current
                
            else:
                
                T_High = T_Current
        Temp_dif = abs(T_High - T_Low)
    
            
    return T_Current                


# In[122]:


def Temp_Range_Equilibrium(c_alloy):
    
    # Function to find Temp range over which Liq + Solid equilibria exists 
    # Liquidus and Solidus Temp for a alloy composition : c_alloy is found
    
    start = 400       # start : Temperature below which only Solid Exists
    stop = 1000       # stop : Temperature above which only Liquid Exists
    
    T_Bin_Res = Binary_Search(start,stop,c_alloy) 
    T_Bin_Res[0:-1] = sorted(T_Bin_Res[0:-1])
    
    if T_Bin_Res[3] == 1: 
        
        T_Liquidus = Binary_Melting(T_Bin_Res[1],T_Bin_Res[2],c_alloy,[Liquid_Phase])
        print('Liquidus Temperature:',T_Liquidus)
        #print(equi_composition_phase(T_Liquidus,c_alloy,[Liquid_Phase]))
        #print(phase_in_equilibrium(T_Liquidus,c_alloy))
        print('')
        T_Solidus =  Binary_Melting(T_Bin_Res[0],T_Bin_Res[1],c_alloy,[Solid_Phase])
        print('Solidus Temperture:',T_Solidus)
        #print(phase_in_equilibrium(T_Liquidus,c_alloy))
        #print(equi_composition_phase(T_Solidus,c_alloy,[Solid_Phase]))
        
        return [T_Solidus,T_Liquidus]
        
    else:
        
        print('Equilibrium Not Found')
        return 


# In[123]:


def Tem_equi_composition(start,stop,c_alloy):
    
    # Function that writes chemical potentials and equilibrium compositions(w.r.t Compo_Calc_1) in the order below:
    # Temp, chem_pot of Compo_Calc_1,chem_pot of Compo_Calc_2, equi_compo(solid), equi_compo(liquid)
    
    # alloy composition: c_alloy
    # start : temperature at which equilibria starts for the alloy composition 
    # stop : temperature at which equilibria ends for the alloy composition
    
    del_ta = 1
    df_Final = df()
    df_der_Solid = df()
    df_der_Liquid = df()
    
    T_array = np.arange(start,stop+del_ta,del_ta) 
    phases_at_equi = sorted(phases_equi)
    
    for Temp in T_array:        
        
        phase_composition_equi = []
        
        eq1 = equilibrium(db_sys, my_component_sys, my_phases_sys, {v.X(Compo_calc_1): c_alloy ,v.T: Temp, v.P: 101325},output='GM')
        Phase_equili_T = np.array(eq1.Phase.squeeze())
        Phase_equili_T = list(Phase_equili_T)
        Phase_equili_T = [x for x in Phase_equili_T if x]
        Phase_equili_T = sorted(Phase_equili_T)
        
        df_Temp = df()
        der_1 = df()   
        der_2 = df()
    
        if phases_at_equi == Phase_equili_T :                      
            
            for phase_equi in phases_equi:    #phases_equi : Ist element is always Solid and last is Liquid
              
                phase_compo_X = np.array(eq1.X.where(eq1.Phase == phase_equi).sel(component = Compo_calc_1).squeeze())
                phase_compo_X = phase_compo_X[~(np.isnan(phase_compo_X))]        
                phase_composition_equi.append(phase_compo_X)
                                
                if phase_equi == Solid_Phase:
                   
                   my_list = [v.Y(phase_equi,0,calc_component_sys[ind_ex_1]), v.Y(phase_equi,0,calc_component_sys[ind_ex_2]), v.T]
                   der_val = d_dG_Sol.subs([(my_list[0], phase_compo_X[0]),(my_list[1], 1-phase_compo_X[0]),(my_list[2], Temp)])
                   der_1 = df({'T':Temp, strt_A: der_val },index=[0])                
                
                if phase_equi == Liquid_Phase:
                
                   my_list = [v.Y(phase_equi,0,calc_component_sys[ind_ex_1]), v.Y(phase_equi,0,calc_component_sys[ind_ex_2]), v.T]
                   der_val = d_dG_Liq.subs([(my_list[0], phase_compo_X[0]),(my_list[1], 1-phase_compo_X[0]),(my_list[2], Temp)])
                   der_2 = df({'T':Temp, strt_B: der_val },index=[0])           
                
            df_Temp = df({'T':Temp,strt_C: phase_composition_equi[0],strt_D: phase_composition_equi[1]},index=[0])
            
            
            df_Final = df_Final.append(df_Temp)
            df_der_Solid = df_der_Solid.append(der_1)
            df_der_Liquid = df_der_Liquid.append(der_2)
                            
    return df_Final, df_der_Solid, df_der_Liquid

# In[]:

def d_dG_eq(phase):
    
    ind_ex_1 = calc_component_sys.index(Compo_calc_1)    
    ind_ex_2 = calc_component_sys.index(Compo_calc_2)  

    print(phase)
    print('')
    phase_mod = 'mod_' + phase
    phase_mod = Model(db_sys, my_component_sys, phase)
    my_list = [v.Y(phase,0,calc_component_sys[ind_ex_1]), v.Y(phase,0,calc_component_sys[ind_ex_2]), v.T]
          
    var = my_list
       
    GM_phase = phase_mod.GM
    GM_phase_diff = {}
    GM_phase_doub_diff = {}
      
    for n in var[0:-1]:
       
        GM_phase_diff[n] = GM_phase.diff(n)
   
    dG_phase = GM_phase_diff[var[ind_ex_1]] - GM_phase_diff[var[ind_ex_2]]
   
    Dou_dg_phase = dG_phase.diff(var[ind_ex_1]) - dG_phase.diff(var[ind_ex_2]) 

    return Dou_dg_phase 
# In[127]:

d_dG_Sol = d_dG_eq(Solid_Phase)
d_dG_Liq = d_dG_eq(Liquid_Phase)

T_Range = Temp_Range_Equilibrium(c_alloy)


# In[128]:


T_Range = sorted(T_Range)
start = T_Range[0]

if (start - int(start)) != 0:
    
    start = int(start) + 1

stop = int(T_Range[1])
#print(start,stop)
df_compo, d2_G_Sol, d2_G_Liq = Tem_equi_composition(start,stop,c_alloy)


# In[130]:


df_compo.to_csv("Equilibrium_Composition.txt", sep = ',',index = False)

d2_G_Sol.to_csv(strt_A, sep = ',',index = False)
d2_G_Liq.to_csv(strt_B, sep = ',',index = False)
# In[133]:


def gen_c_phase_eq(phase_equil):
     
   # Function to generate c codes with Free Energy(G), derivative G and double derivative of G w.r.t compo_calc_1     
        
    print('Generating .c files')       
    print('')
    
    for phase in phase_equil:        
        
        print(phase)
        print('')
        phase_mod = 'mod_' + phase
        phase_mod = Model(db_sys, my_component_sys, phase)
                        
        y = sym.MatrixSymbol('y', phase_mod.GM.args.__len__()-1, 1)
        T = sym.MatrixSymbol('T', 1, 1)
        var = tuple(y)+tuple(T)
        my_list = [v.Y(phase,0,calc_component_sys[ind_ex_1]), v.Y(phase,0,calc_component_sys[ind_ex_2]), 'T']
               
        mapped = zip(my_list, list(var))
        mapped = list(mapped)
        state_array_map= dict(zip(my_list, list(var)))
     #   print(state_array_map)
            
        GM_phase = 'reduced_GM_' + phase           
        
        GM_phase = phase_mod.GM.xreplace(state_array_map)
        
        GM_phase_diff = 'reduced_GM_diff_' + phase
        
        GM_phase_diff = {}
        GM_phase_doub_diff = {}
           
        for n in var[0:-1]:
            
            GM_phase_diff[n] = GM_phase.diff(n)
        
        dG_phase = GM_phase_diff[var[ind_ex_1]] - GM_phase_diff[var[ind_ex_2]]
        
        Dou_dg_phase = dG_phase.diff(var[ind_ex_1]) - dG_phase.diff(var[ind_ex_2]) 
        
        G = sym.Symbol('dG')
        dG = sym.Symbol('dG')
        d2_G = sym.Symbol('d2_G')
        eqn_G = sym.Eq(G, GM_phase)
        eqn = sym.Eq(dG, dG_phase)
        eqn1 = sym.Eq(d2_G, Dou_dg_phase)
        
        if phase == Liquid_Phase:
          
            strt_ther = 'G_' +  phase 
            strt_thermo = 'DG_' +  phase  
            strt_chempot = 'dou_DG_' + phase
    
        else:
          
            strt_ther = 'G_SOLID' 
            strt_thermo = 'DG_SOLID'  
            strt_chempot = 'dou_DG_SOLID' 
    
        codegen((strt_ther, eqn_G), language='c', to_files=True)
        
        codegen((strt_thermo, eqn), language='c', to_files=True)

        codegen((strt_chempot, eqn1), language='c', to_files=True)


# In[134]:


from sympy.utilities.codegen import codegen
gen_c_phase_eq(phases_equi)


# In[ ]:





# In[ ]:




