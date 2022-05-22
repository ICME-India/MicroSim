#!/bin/bash

echo "Start"


python FE_Equi_Composition_Original.py alzn_mey.tdb FCC_A1 LIQUID ZN AL 0.7 

gcc -c G_SOLID.c
gcc -c G_LIQUID.c
gcc -c DG_SOLID.c
gcc -c DG_LIQUID.c
gcc -c dou_DG_SOLID.c
gcc -c dou_DG_LIQUID.c

echo "Done "



