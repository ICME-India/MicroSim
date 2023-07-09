#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 21:46:00 2021

@author: abhik
"""

## Install pycalphad and related dependencies. Sympy gets automatically installed
## You will also need pandas

import string
from tinydb import where
import sympy as sym

sym.init_printing()

import matplotlib.pyplot as plt
from pycalphad import Model, Database, calculate, equilibrium, binplot
import pycalphad.variables as v

import pandas as pd
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import UnivariateSpline

from scipy.interpolate import CubicSpline

import numpy as np

from sympy import *

from sympy.codegen.ast import Assignment

import sys
import re
import os 

## Generates thermo files for CPU
print('Reading Input file')
with open(sys.argv[1], 'r') as my_file:
    flines=my_file.readlines()
    flag = 0
    index = 0
    for line in flines:
        index += 1
        searchstatus = re.match(r'\btdbfname\b', line)
        if searchstatus:
            flag = 1
            break
    if flag == 0:
        print("###########################################")
        print("# Error in Input file: No TDB information #")
        print("###########################################")
    else:
        tdbfline = flines[index-1]
        tdbfname = tdbfline.split("=")[1].replace(";", "").strip()
        #print(tdbfname)

    flag = 0
    index = 0
    for line in flines:
        index += 1
        searchstatus = re.match(r'\bCOMPONENTS\b', line)
        if searchstatus:
            flag = 1
            break
    if flag == 0:
        print("##################################################")
        print("# Error in Input file: No components information #")
        print("##################################################")
    else:
        COMPONENTSline = flines[index-1]
        components1 = COMPONENTSline.split("=")[1].replace(';', '').strip()
        components = re.sub(r"[{} ]","",components1).split(",")
        #print(components)

    flag = 0
    index = 0
    for line in flines:
        index += 1
        searchstatus = re.match(r'\bPHASES\b', line)
        if searchstatus:
            flag = 1
            break
    if flag == 0:
        print("#############################################")
        print("# Error in Input file: No phase information #")
        print("#############################################")
    else:
        nPHASESline = flines[index-1]
        nphases1 = nPHASESline.split("=")[1].replace(';', '').strip()
        nphases = re.sub(r"[{} ]","",nphases1).split(",")
        
    flag = 0
    index = 0
    for line in flines:
        index += 1
        searchstatus = re.match(r'\btdb_phases\b', line)
        if searchstatus:
            flag = 1
            break
    if flag == 0:
        print("#############################################")
        print("# Error in Input file: No phase information #")
        print("#############################################")
    else:
        PHASESline = flines[index-1]
        phases1 = PHASESline.split("=")[1].replace(';', '').strip()
        phases = re.sub(r"[{} ]","",phases1).split(",")
        #print(phases)
        
comps=components+['VA']
#print(comps)

if 'fcc'==phases[0].lower():
    phases[0]='FCC_A1'
if 'bcc'==phases[0].lower():
    phases[0]='BCC_A2'
if 'hcp'==phases[0].lower():
    phases[0]='HCP_A3'
if 'liquid' == phases[1].lower() or 'liq' == phases[1].lower():
    phases[1]='LIQUID'

print("\n````````````````````````````````")
print("TDB file name =", tdbfname)
print("Components =", comps)
print("tdb_phases =", phases)
print("````````````````````````````````\n")


##############################  Writing files for OpenCL #defines

nphadef = '#define npha (' + str(len(nphases)) + ')\n'
ncomdef = '#define ncom (' + str(len(components)) + ')\n'
nsoldef = '#define nsol (ncom-1)\n'
#nsoldef = '#define nsol (' + str(len(components-1)) + ')\n'

print(nphadef)
print(ncomdef)

file = open("solverloop/defines.h", "w")
file.write(nphadef)
file.write(ncomdef)
file.write(nsoldef)

# with fileinput.FileInput('solverloop/defines.h', inplace=True) as file:
#     for line in file:
#         if line.strip().startswith('#define npha ('):
#             line = nphadef
#         sys.stdout.write(line)


# with fileinput.FileInput('solverloop/defines.h', inplace=True) as file:
#     for line in file:
#         if line.strip().startswith('#define ncom ('):
#             line = ncomdef
#         sys.stdout.write(line)


###############################################################


#my_phases_tdb = ['LIQUID', 'FCC_A1']
### Reading in the .tdb files
tdbf = Database(tdbfname)

## You can initialize a phase string consisting of the phases. The list you can
## by just printing tdbf.phases.keys(). You can just create a list of the phases
## that you want to consider. If you want all the phases you can also do
## my_phases_tdb = list(tdbf.phases.keys())

#my_phases_tdb = ['LIQUID', 'FCC_A1']
#comps = ['Nb', 'Ni', 'VA']
#print(len(comps))
## Importing the codegen library for codeprinting
from sympy.utilities.codegen import codegen

def generate_thermodynamic_files(tdbf, comps, phase):
  ## Initializing the instance of the model for a particular phase in this
  ## case 'FCC_A1', It is has different attributes, you can just do
  ## dir(mod3) to print out the possibilities. Here we will utilize
  ## the free-energies that are given by mod3.GM
  mod3 = Model(tdbf, comps, phase)
  print("\nGE =", mod3.GM,"\n")
  
  ##You will see that the functions have long variable names
  ## corresponding to the compositions of the sublattice compositions.
  ## In order to create a clean C-code, we will replace the variables
  ## by performing a map between the new variables and the old variables

  ## These are the set of new variables that are being created
  ## This is held in 'var' which is essentially the symbols for the
  ## two sublattices concentrations of 'comps' and followed by T
  y = sym.MatrixSymbol('y', mod3.GM.args.__len__()-1, 1)
  T = sym.MatrixSymbol('T', 1, 1)
  var = tuple(y)+tuple(T)
  print(mod3.GM.args.__len__())
  ## The original set of variables in mod3.GM. The sublattice
  ## variables are obtained from v.Y(phase_name, sublattice_index, species_name)
  ## Derivatives are calculated w.r.t. arguments in below lists excluding temperature
  #my_list = [v.Y(phase,0,comps[0]), v.Y(phase,0,comps[1]), 'T']
  #sym.symbols('x_%d' % i) for i in range(3)
  #my_list = [v.Y('FCC_A1',0,'AL'), v.Y('FCC_A1',0,'ZN'), 'T']
  my_list = [v.Y(phase, 0, comps[i]) for i in range(len(components))]
  my_list += T.name
  print("my_List -> ", my_list)
  ## Creating a mapping between the old and the new
  ## that is stored in state_array_map

  mapped = zip(my_list, list(var))
  mapped = list(mapped)
  state_array_map= dict(zip(my_list, list(var)))
  print("Mapped -> ", state_array_map)
  
  # Replacing the original variables in the mod3.GM expression
  # with the state_array_map. The new expression is stored in
  # a new expression reduced_GM

  reduced_GM = mod3.GM.xreplace(state_array_map)
  print("\nGE_%s =" % phases[0], reduced_GM,"\n")
  
  # Deciding the 'var' boundaries for partial derivatives
  # Presence of magnetic contribution may alter the size of 'var'

  if (len(var)==len(comps)):
      a1=1
  elif ((len(var)-1)==len(comps)):
      #magnetic contribution might be present
      a1=2
  else:
      print("Unknown vribales in tuple var")
      print("Expressions may not be correct")
      #print("Manual editing of this file is necessary")
      print("No expressions are generated")
      print("Presently, code works for only binary system")
      print("Exiting")
      exit()
  print("a1=",a1)
  
  ## Doing the same for the partial derivatives of the
  ## free-energy vs sublattice compositions.
  ## The partial derivatives are being stored here in a
  ## dictionary with the names of the new variables as keys.
  ## Derivatives are stored in reduced_GM_diff
  ## Note, that derivatives can simply be taken by
  ## .diff operation

  print('Calculating partial derivatives')
  reduced_GM_diff={}
  for i, n in enumerate(var[0:len(components)-1]):
      reduced_GM_diff[n] = reduced_GM.diff(n) - reduced_GM.diff(var[len(components)-1])

  #Converting the dictionary of the derivatives into a list

  mu = [reduced_GM_diff[n] for n in var[0:len(components)-1]]

  # And into a matrix. This is required for codegen
  # for printing a matrix eqn into a C-code

  mu_expr = sym.Matrix(mu)
  print(mu_expr)
  
  print('Calculating second derivatives')
  reduced_GM_diff2={}
  for i, n1 in enumerate(var[0:len(components)-1]):
    reduced_GM_diff2[n1]={}
    
  for i, n1 in enumerate(var[0:len(components)-1]):
    for j, n2 in enumerate(var[0:len(components)-1]):
      reduced_GM_diff2[n1][n2] = mu[i].diff(n2) - mu[i].diff(var[len(components)-1])

  #Converting the dictionary of the derivatives into a list

  dmudc = [[reduced_GM_diff2[n1][n2] for n2 in var[0:len(components)-1]] for n1 in var[0:len(components)-1]]

  # And into a matrix. This is required for codegen
  # for printing a matrix eqn into a C-code

  dmudc_expr = sym.Matrix(dmudc)
  
  #print(dmudc_expr[0,0]);
  
  return reduced_GM, mu_expr, dmudc_expr, len(components), len(phases)


free_energy_tdb = [sym.symbols('GE_%d' %i)    for i in range(len(phases))]
Mu_tdb          = [sym.symbols('Mu_%d' %i)    for i in range(len(phases))]
dmudc_tdb       = [sym.symbols('dmudc_%d' %i) for i in range(len(phases))]


eqn_GM          = [sym.symbols('eqn_GM_%d' %i)    for i in range(len(phases))]
eqn_mu          = [sym.symbols('eqn_mu_%d' %i)    for i in range(len(phases))]
eqn_dmudc       = [sym.symbols('eqn_dmudc_%d' %i) for i in range(len(phases))]


for i in range(len(phases)):
  GM, mu, dmu_dc, num_comps, num_tdb_phases = generate_thermodynamic_files(tdbf, comps, phases[i])
  
  GE           = sym.symbols('Ge')
  MU           = sym.MatrixSymbol('Mu',    num_comps-1,1)
  DMUDC        = sym.MatrixSymbol('Dmudc', num_comps-1,num_comps-1)
  
  eqn_GM[i]    = sym.Eq(GE,    GM)
  eqn_mu[i]    = sym.Eq(MU,    mu)
  eqn_dmudc[i] = sym.Eq(DMUDC, dmu_dc)
  

functions  = [(free_energy_tdb[i].name, eqn_GM[i]) for i in range(len(phases))]
functions += [(Mu_tdb[i].name, eqn_mu[i])          for i in range(len(phases))]
functions += [(dmudc_tdb[i].name, eqn_dmudc[i])    for i in range(len(phases))]

print(functions)

if not os.path.exists('tdbs'):
    os.makedirs('tdbs')
#os.mkdir('tdbs', mode=0o777)

codegen(functions, language='c', to_files=True, prefix='tdbs/Thermo')

import fileinput

with fileinput.FileInput('tdbs/Thermo.c', inplace=True) as file:
    for line in file:
        print(line.replace('#include <math.h>', '//#include <math.h>'), end='')


free_energy_functions = "{" + free_energy_tdb[0].name

for i in range(1,len(phases)):
 free_energy_functions += "," + free_energy_tdb[i].name

free_energy_functions += "}"

Diffusion_potential_functions = "{" + Mu_tdb[0].name

for i in range(1,len(phases)):
 Diffusion_potential_functions += "," + Mu_tdb[i].name

Diffusion_potential_functions += "}"

Hessian_functions = "{" + dmudc_tdb[0].name 

for i in range(1,len(phases)):
 Hessian_functions += "," + dmudc_tdb[i].name

Hessian_functions += "}"

#print(str);
## For .h File: void... added before #endif
with fileinput.FileInput('tdbs/Thermo.h', inplace=True) as file:
    for line in file:
        print(line.replace('#endif', 'void(*free_energy_tdb[])(double T, double *y, double *Ge) = %s; \n#endif '% free_energy_functions), end='')

with fileinput.FileInput('tdbs/Thermo.h', inplace=True) as file:
    for line in file:
        print(line.replace('#endif', 'void(*Mu_tdb[])(double T, double *y, double *Mu) = %s; \n#endif '% Diffusion_potential_functions), end='')

with fileinput.FileInput('tdbs/Thermo.h', inplace=True) as file:
    for line in file:
        print(line.replace('#endif', 'void(*dmudc_tdb[])(double T, double *y, double *Dmudc) = %s; \n#endif '% Hessian_functions), end='')


# Generating Files for OpenCL 



codegen(functions, language='c', to_files=True, prefix='solverloop/ThermoCL_2')

import fileinput

with fileinput.FileInput('solverloop/ThermoCL_2.c', inplace=True) as file:
    for line in file:
        print(line.replace('#include <math.h>', '//#include <math.h>'), end='')



#################### Writing thermo function calls for OpenCL

file = open("solverloop/GibbsEnergyData2.h", "w")
file.write("//Do not edit this file alone. Function calls in kerenls can go wrong if any changes are made...\n\n")
file.write("#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n")
file.write("void Ge(double T, double *y, double *Ge, int tpha);\n")
file.write("void Mu(double T, double *y, double *MU, int tpha);\n")
file.write("void dMudc(double T, double *y, double *Dmudc, int tpha);\n")
file.write("\n\n")
file.write("void Ge(double T, double *y, double *Ge, int tpha) { \n")

Ge_functions = "\n  if ( tpha == 0 ) { \n " + "   " + free_energy_tdb[0].name + "(T, y, Ge); \n  }\n"

for i in range(1,len(phases)):
    Ge_functions += "  else if ( tpha == " + str(i) +" ) { \n"
    Ge_functions += "    " + free_energy_tdb[i].name + "(T, y, Ge); \n  }\n"

Ge_functions += "\n}\n\n"

file.write(Ge_functions)

file.write("void Mu(double T, double *y, double *MU, int tpha) { \n")
Mu_functions = "\n  if ( tpha == 0 ) { \n " + "   " + Mu_tdb[0].name + "(T, y, MU); \n  }\n"

for i in range(1,len(phases)):
    Mu_functions += "  else if ( tpha == " + str(i) +" ) { \n"
    Mu_functions += "    " + Mu_tdb[i].name + "(T, y, MU); \n  }\n"

Mu_functions += "\n}\n\n"

file.write(Mu_functions)

file.write("void dMudc(double T, double *y, double *Dmudc, int tpha) { \n")
dMudc_functions = "\n  if ( tpha == 0 ) { \n " + "   " + dmudc_tdb[0].name + "(T, y, Dmudc); \n  }\n"

for i in range(1,len(phases)):
    dMudc_functions += "  else if ( tpha == " + str(i) +" ) { \n"
    dMudc_functions += "    " + dmudc_tdb[i].name + "(T, y, Dmudc); \n  }\n"

dMudc_functions += "\n}\n\n"

file.write(dMudc_functions)
