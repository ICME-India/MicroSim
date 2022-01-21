#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 21:46:00 2021

@author: abhik
"""
## This code generates GLiq.h GLiq.c GSol.h and GSol.c files at solverloop/GibbsEnergyFunctions/
## Above files are necessesary for getting Gibbs energy and their derivative values in 
## Kim solidification module 

## Install pycalphad and related dependencies. Sympy gets automatically installed
## You will also need pandas

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
if 'fcc'==phases[1].lower():
    phases[1]='FCC_A1'
if 'bcc'==phases[1].lower():
    phases[1]='BCC_A2'
if 'hcp'==phases[1].lower():
    phases[1]='HCP_A3'

print("\n````````````````````````````````")
print("TDB file name =", tdbfname)
print("Components =", comps)
print("Phases =", phases)
print("````````````````````````````````\n")

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

## Initializing the instance of the model for a particular phase in this 
## case 'FCC_A1', It is has different attributes, you can just do 
## dir(mod3) to print out the possibilities. Here we will utilize 
## the free-energies that are given by mod3.GM
mod3S = Model(tdbf, comps, phases[0])
print("\nGE =", mod3S.GM,"\n")

##You will see that the functions have long variable names 
## corresponding to the compositions of the sublattice compositions.
## In order to create a clean C-code, we will replace the variables 
## by performing a map between the new variables and the old variables

## These are the set of new variables that are being created
## This is held in 'var' which is essentially the symbols for the 
## two sublattices concentrations of 'comps' and followed by T
y = sym.MatrixSymbol('y', mod3S.GM.args.__len__()-1, 1)
T = sym.MatrixSymbol('T', 1, 1)
var = tuple(y)+tuple(T)
#print(len(var))

## The original set of variables in mod3.GM. The sublattice 
## variables are obtained from v.Y(phase_name, sublattice_index, species_name)
## Derivatives are calculated w.r.t. arguments in below lists excluding temperature
my_listS = [v.Y(phases[0],0,comps[0]), v.Y(phases[0],0,comps[1]), 'T']
my_listL = [v.Y(phases[1],0,comps[0]), v.Y(phases[1],0,comps[1]), 'T']

print("my_ListS -> ", my_listS)

## Creating a mapping between the old and the new
## that is stored in state_array_map
mapped = zip(my_listS, list(var))
mapped = list(mapped)
state_array_map= dict(zip(my_listS, list(var)))
print("Mapped -> ", state_array_map)

## Replacing the original variables in the mod3.GM expression 
## with the state_array_map. The new expression is stored in 
## a new expression reduced_GM
reduced_GMS = mod3S.GM.xreplace(state_array_map)
print("\nGE_%s =" % phases[0], reduced_GMS,"\n")

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

#print("a1=",a1)

### Doing the same for the partial derivatives of the 
### free-energy vs sublattice compositions.
### The partial derivatives are being stored here in a
### dictionary with the names of the new variables as keys.
### Derivatives are stored in reduced_GM_diff
### Note, that derivatives can simply be taken by 
### .diff operation
print('Calculating partial derivatives')
reduced_GMS_diff={}
for i, n in enumerate(var[0:-a1]):
    reduced_GMS_diff[n] = reduced_GMS.diff(n)

##Converting the dictionary of the derivatives into a list
expr1 = [reduced_GMS_diff[n] for n in var[0:-a1]]

## And into a matrix. This is required for codegen 
## for printing a matrix eqn into a C-code
mat_expr = sym.Matrix(expr1)

## We are now going to be creating equations which can be written 
## into C-code, as functions. The lhs of the equations are going to 
## symbols or Matrix symbols depending on whether you want to write 
## single equation or multiple equations

## This symbol in short for Gibbs energy of solid, will be used to denote
## the lhs of the Gibbs energy expression
GES = sym.symbols('GES')

## This Matrix symbol (dGES) is used for storing the lhs for the 
## partial derivatives of the Gibbs energies with composition
## Note the shape of the matrix is the same as the number of 
## partial derivatives
dGES = sym.MatrixSymbol('dGES', mod3S.GM.args.__len__()-a1,1)

## Creating the equations using the symbolic equation 
## creator sym.Eq
eqn_GES = sym.Eq(GES, reduced_GMS)
eqn_dGES = sym.Eq(dGES, mat_expr)

muS = expr1[0] - expr1[1]

reduced_GMS_diff2={}
for i, n in enumerate(var[0:-a1]):
    reduced_GMS_diff2[n] = muS.diff(n)


expr1a = [reduced_GMS_diff2[n] for n in var[0:-a1]]

mat_expr1 = sym.Matrix(expr1a)

ddGES = sym.MatrixSymbol('ddGES', mod3S.GM.args.__len__()-a1,1)

eqn_ddGES = sym.Eq(ddGES, mat_expr1)

## Writing the equations to the files specified by 'prefix'

## This creates GSol.c and GSol.h, where GSol.c
## contains the functions of Gibbs energy and partial derivatives, whose output is returned as 
## a pointer that is passed as an argument to the function
codegen([('GES', eqn_GES), ('dGES', eqn_dGES), ('ddGES', eqn_ddGES)], language='c', to_files=True, prefix='solverloop/GibbsEnergyFunctions/GSol')

###############################################
##Phase 2 descriptions
print('\n##',phases[1],'Gibbs energy ##\n')
mod3L = Model(tdbf, comps, phases[1])
#print(mod3L.GM)

y = sym.MatrixSymbol('y', mod3L.GM.args.__len__()-1, 1)
T = sym.MatrixSymbol('T', 1, 1)
var = tuple(y)+tuple(T)
print("my_ListL -> ", my_listL)

mapped = zip(my_listL, list(var))
mapped = list(mapped)
state_array_map= dict(zip(my_listL, list(var)))
print("Mapped -> ", state_array_map)

reduced_GML = mod3L.GM.xreplace(state_array_map)
print("\nGE_%s =" % phases[1], reduced_GML,"\n")

print('Calculating partial derivatives')
a1l=1
reduced_GML_diff={}
for i, n in enumerate(var[0:-a1l]):
    reduced_GML_diff[n] = reduced_GML.diff(n)


expr1L = [reduced_GML_diff[n] for n in var[0:-a1l]]

mat_exprL = sym.Matrix(expr1L)

GEL = sym.symbols('GEL')

dGEL = sym.MatrixSymbol('dGEL', mod3L.GM.args.__len__()-a1l,1)

eqn_GEL = sym.Eq(GEL, reduced_GML)
eqn_dGEL = sym.Eq(dGEL, mat_exprL)

muL = expr1L[0] - expr1L[1]

reduced_GML_diff2={}
for i, n in enumerate(var[0:-a1l]):
    reduced_GML_diff2[n] = muL.diff(n)

expr1aL = [reduced_GML_diff2[n] for n in var[0:-a1l]]

mat_expr1L = sym.Matrix(expr1aL)

ddGEL = sym.MatrixSymbol('ddGEL', mod3L.GM.args.__len__()-a1l,1)

eqn_ddGEL = sym.Eq(ddGEL, mat_expr1L)

codegen([('GEL', eqn_GEL), ('dGEL', eqn_dGEL), ('ddGEL', eqn_ddGEL)], language='c', to_files=True, prefix='solverloop/GibbsEnergyFunctions/GSol_m')

##########################
import fileinput

with fileinput.FileInput('solverloop/GibbsEnergyFunctions/GSol_m.c', inplace=True) as file:
    for line in file:
        print(line.replace('#include <math.h>', '//#include <math.h>'), end='')


with fileinput.FileInput('solverloop/GibbsEnergyFunctions/GSol.c', inplace=True) as file:
    for line in file:
        print(line.replace('#include <math.h>', '//#include <math.h>'), end='')

print("\nFiles generated successfully\n")
