#!/bin/bash

if [ $# -ne 3 ]; 
then 
	echo "Arguments are missing"
	echo "Make sure arguments are Input file, Filling file and Output name"
	echo "For example arguments can be"
	echo "Input.in Filling.in Output"
	exit
else
	python3 GEdata_writer.py $1

	./kim_soldfn.out $1 $2 $3
fi
