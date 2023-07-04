#!/bin/bash

make

g++ -o Replace Replace.cpp

./Replace Input.in Filling.in
#./Replace Input_TP.in Filling_TP.in

mpirun -np 4 ./main2d.gnu.MPI.ex input2.in
