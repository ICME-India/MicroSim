# MicroSim
MicroSim is a software stack that consists of phase-field codes that offer flexibility with discretization, models as 
well as the high-performance computing hardware(CPU/GPU) that they can execute on. Along with this the stack also consists of Multi-physics solver modules 
that are based on OpenFoam and AMRex libraries(will be added soon). The stack has an integrator interface that is built 
using python that allows one to create the input and filling files required for the solvers as well as provides 
a consolidated framework to choose the solver, compile, execute and visualize simulation results.
The project is a consortium between (IISc Bangalore, IIT Hyderabad, IIT Bombay, IIT Madras, Savitribai Phule Pune University, C-DAC Pune).
Following is a brief description of the different software modules and details on how to independently execute them.  

 [Workshop-cum-demo](https://www.youtube.com/channel/UCnmLKb-kqQCNXWa0Oz6Ua_g)  
 
## Grand-Potential Model(SERIAL)
This is multiphase multi-component phase-field solver based on the 
grand-potential formalism. This is a 2D serial CPU version of the code.


Compilation can be done by a simple "make".


For running the following execution command is required

./microsim_gp name_of_infile name_of_filling_file name_of_output_files

The files are written in the DATA folder in the .vtk format.

The code has been developed at IISc Bangalore, Department of Materials 
Engineering by Prof. Abhik Choudhury. 

## Grand-Potential Model(MPI)

This is multiphase multi-component phase-field solver based on the 
grand-potential formalism. The solver is parallelized using MPI on 
CPUs and requires the h5pcc compiler for compilation and execution. 

Compilation can be done by a simple "make"

For running the code on the cluster please use a script. 
For testing on your desktops/laptops the following execution command is required

mpirun -np num_processors ./microsim_gp name_of_infile name_of_filling_file name_of_output_files num_workers_x num_workers_y

For .h5 files, with WRITEHDF5=1, output files need to be transformed in .xml format using the following command
just above the DATA folder that is created upon execution

./write_xdmf name_of_infile name_of_output_file total_no_workers start_time end_time

For ASCII files in .vtk format the consolidated output files needs to be reconstructed out of separate processor files
that are written in the DATA folder that is created upon execution

./reconstruct name_of_infile name_of_output_file number_of_workers start_time end_time


The code has been developed at IISc Bangalore, Department of Materials 
Engineering by Prof. Abhik Choudhury. 

The code is built on work by past PhD students 

a) Sumeet Rajesh Khanna


## Cahn–Hilliard Model
 * C code for precipitate evolution. 
 * It solves Allen-Cahn and Cahn–Hilliard equations using FFTW3. 
 * 
 * Compile the code using "make". Compilation creates "FFT_2D_ppt.out" file. 
 * Execute "./microsim_ch_fft Input.in Filling.in output" 
 * 
 * Authors: Dasari Mohan and M P Gururajan
 * 
 * This is alpha version of the code and check for updates in future release. 
 * 
 * Copyright (c) 2021 Materials and Process Modelling Laboratory,
 * Department of Metallurgical Engineering and Materials Science, 
 * Indian Institute of Technology Bombay, Mumbai 400076 INDIA.

## KKS GPU CUDA Model
This module solves a multiphase, multicomponent KKS-formulation of the phase-field equations.
The code is written using CUDA and OpenMPI, and can be employ multiple GPUs in a multi-node setting.
It has been tested using Nvidia Tesla V100s, Nvidia Tesla P100s, and Nvidia Tesla K80s, with CUDA 11.x, 12.x and OpenMPI 4.0.x.

Follow the instructions in this module's README to get this solver up and running.

To compile the solver, simply open a terminal in the base directory of the module and run the command 'make'.
$ make
To compile without cuFFTMp or HDF5 support,
$ make ENABLE_CUFFTMP=0 ENABLE_HDF5=0

For usage on the PARAM supercomputers, one can use the SLURM script (ParamJobScript.sh) and Makefile (Makefile_Param) that are included.
Since the packages may differ from platform to platform, some modifications to the above may be necessary.

To run the solver, use:
$ make run INPUT=<name_of_infile> FILLING=<name_of_filling_file> OUTPUT=<name_of_output_file> NPROCS=<number_of_processors>
or

 - GPU Phase-Field Developer Team @ IITH (Saurav Shenoy, Saswata Bhattacharya)

The following contributors are acknowledged
    Tushar Jogi
    Pankaj
    Hemanth Kumar Sandireddy


## KKS GPU OPENCL Model
 * OpenCL code for solidification microstructure evolution 
 * 
 * Compile the code using "make". Compilation creates "kim_soldfn.out" file. 
 * GEdata_writer.py is used for generation of Gibbs energy and its derivatives
 * To generate Gibbs energies and execute the program
 * run "./kimsldfn.sh  Input.in Filling.in Output"
 * It is always safe to run above command for execution of the code.
 *
 * If Gibbs energies are generated already then generating 
 * Gibbs energies can be skipped and directly execute following command.
 * Execute "./microsim_kks_opencl Input.in Filling.in Output" 
 * 
 * Authors: Dasari Mohan and G Phanikumar
 * Acknowledgement to P. Gerald Tennyson for contributions towards code development at IITM
 * 
 * This is alpha version of the code and check for updates in future release. 
 * 

## Grand-potential OpenFOAM

An OpenFOAM phase-field solver to simulate solidification of binary and ternary alloys

### solver

It contains the source files of the solver.

### cases
#### multigrainAlZn

It contains the OpenFOAM case files required to run the multigrain problem for AlZn alloy.

#### coolingAlZn

It contains the OpenFOAM case files required to run the cooling problem for AlZn alloy.

#### multigrainNiNb

It contains the OpenFOAM case files required to run the multigrain problem for NiNb alloy.

#### coarseningAlZn

It contains the OpenFOAM case files required to run the coarsening problem for AlZn alloy.

#### multigrainNiAlMo

It contains the OpenFOAM case files required to run the multigrain problem for NiAlMo alloy.

#### coolingNiAlMo

It contains the OpenFOAM case files required to run the cooling problem for NiAlMo alloy.

#### coarseningNiAlMo

It contains the OpenFOAM case files required to run the coarsening problem for NiAlMo alloy.

#### coolingCoarseningNiAlMo

It contains the OpenFOAM case files required to run the coarsening problem while cooling for NiAlMo alloy.


The following contributers are acknowledged
1. Swapnil Bhure
2. Tanmay Dutta 
3. Ravi Kumar Singh
4. Bhalchandra Bhadak


## Infile Generator
Python GUI application for generating Infile and Filling files.

* This script depends on gnome-terminal. So, make sure you have it installed from the package repository of your linux distribution.

* You can use a package manager like Miniconda or Anaconda to avoid issues with system python. Miniconda is enough for this specific purpose.

* Install Miniconda package manager from https://repo.anaconda.com/miniconda/Miniconda3-py310_23.1.0-1-Linux-x86_64.sh

* Now create a virtual environment with python 3.9 (version previous to this are also compatible with the packages required) and pip:

> conda create --name msenv python=3.9 pip

* Activate virtual environment msenv:

> conda activate msenv

* Install the packages below using pip (group them together to avoid dependency issue):

> pip install pyqt5 scikit-image vtk tinydb sympy==1.8 pycalphad==0.9.2 pymks yt

* Launch MicroSim:

> python MicroSim.py

### Do you want to modify the GUI as a developer?

* You can switch to some other environment like base:

> conda activate base

* Install pyqt to have access to QT designer:

> conda install pyqt

* Launch QT designer:

> designer

* Open resources/mainscreen.ui
     
Developed by- Ajay Sagar and Tanmay Dutta
