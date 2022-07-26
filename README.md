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
This code solves the problem of precipitate growth using a multiphase field solver on the GPU.
The code is tested on Tesla P100 and Tesla V100.
For Tesla K80, one needs to comment "CudaMemPrefetchAsync" in solverloop/evolve.h

To run the code use 
1. "make" to create executable microsim_kks_cufft
2. Execute "./microsim_kks_cufft Input.in Filling.in Output"

The code uses CUDA version 11, CUFFT and CUDA-CUB libraries. 
nvcc version 11.2 is used for compilation of the codes.
Input.in contains all numerical and physical parameters used in the simulations. 
Makefile creates a  DATA folder where the datafiles (in the form of VTK files) are stored.
VTK files can be viewed using Paraview.

This is the alpha version of code. We will continue to add more features in future release.


- GPU Phase-Field Developer Team @ IITH
  (Pankaj, Saurav Shenoy, Saswata Bhattacharya)

The following contributers are acknowledged
1. Tushar Jogi
2. Hemanth Kumar Sandireddy

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

## MultiPhysics Solver

# Requirements

To build and run the solver modules, the users need OpenFOAM and GSL. The modules are tested using OpenFOAM-in-Box-20 in Ubuntu 18.04 and Ubuntu 20.04, and OpenFOAM-6 in Ubuntu 18.04. Installation image of Ubuntu 20.04 can be downloaded from the link below:

* https://releases.ubuntu.com/20.04.4/ 

OpenFOAM-in-Box-20 is recommended for new users due to ease of installation, and can be downloaded from the links below:

* Official link with installation guide: https://www.cfdsupport.com/openfoam-in-box.html

* Unofficial link: https://drive.google.com/file/d/17gkbpQTK54Hq1A7GNk-s6rktQDfrtcnn/view?usp=sharing

The following commands can be used in the terminal to begin the installation (check official link above if needed):

> mkdir -p ~/OpenFOAM/OpenFOAM-in-Box

> cp OpenFOAM-in-Box-20.09v2-linux64.sh ~/OpenFOAM/OpenFOAM-in-Box/

> cd ~/OpenFOAM/OpenFOAM-in-Box/

> bash OpenFOAM-in-Box-20.09v2-linux64.sh -install

> echo "source ~/OpenFOAM/OpenFOAM-in-Box/OpenFOAM-in-Box-20.09v2-22-g178c07ee/OpenFOAM-dev/etc/bashrc" >> ~/.bashrc

> source ~/.bashrc

Cases can be run after creating the directory below:

> mkdir -p $FOAM_RUN

GSL can be installed using the following command in the terminal:

> sudo apt install libgsl-dev

For visualization and post-processing, ParaView 5.6 and 5.8 are tested. ParaView 5.8 comes along with OpenFOAM in Box 20 and can be launched using the command:

> paraview

To view the plots, gnuplot can be used, which comes preinstalled with Ubuntu. It can be launched using the command:

> gnuplot

It is advised to refer to the OpenFOAM documentions to understand the methods involved while using. The below link can be useful:

* https://cfd.direct/openfoam/documentation/


# OpenFOAM modules

Copy the modules to the OpenFOAM run directory. For instance, to copy the directory PFBinary:

> cp -r PFBinary $FOAM_RUN

Solver has to be compiled from the corresponding solver directory. For instance:

> cd $FOAM_RUN/PFBinary/PFOFBinaryThermo

> wclean

> wmake

Finally, switch to the case directory, e.g. binaryCooling:

> cd $FOAM_RUN/PFBinary/binaryCooling

To run the case with the default parameters, execute the following:

> ./Allclean

> ./Allrun

Note: Allclean and Allrun must be set as executables


## Infile Generator
Python GUI application for generating Infile and Filling files.

### System Requirements :-
   1) Python v3.5 and later 
   2) pyqt5 Library


### Check the Python version on Windows, Linux, or macOS systems :-
To check the version installed, open a terminal window and entering the following:
   
     python --version
![image](https://user-images.githubusercontent.com/20612971/136605964-7aa1de15-1474-4e2f-963b-c3a28c3c6cb1.png)

If installed python version is 3.5 or above then you can continue to install pyqt5 library using pip.
To install python or upgrade older python version follow next step -

### Install Python3

#### Windows
Download setup from https://www.python.org/downloads/ and follow the instruction.
#### linux/Ubuntu
If you are using Ubuntu 16.10 or newer, then you can easily install Python 3 with the following commands:

     sudo apt-get update
     sudo apt-get install python3
If you’re using another version of Ubuntu (e.g. the latest LTS release) or you want to use a more current Python, we recommend using the deadsnakes PPA to install Python 3:

    sudo apt-get install software-properties-common
    sudo add-apt-repository ppa:deadsnakes/ppa
    sudo apt-get update
    sudo apt-get install python3
If you are using other Linux distribution, chances are you already have Python 3 pre-installed as well. If not, use your distribution’s package manager. For example on Fedora, you would use dnf:

    sudo dnf install python3
   
### Download and Install pip:
#### Windows
Download the https://bootstrap.pypa.io/get-pip.py file and store it in the same directory as python is installed.
Run the command given below:
     python get-pip.py
and wait through the installation process.

#### Linux/Ubuntu
To Install PIP3(For Python 3+):

     sudo apt install python3-pip 

### Upgrading Pip version in linux
     sudo -H pip3 install --upgrade pip


### Installing pyqt5 using pip
     pip install pyqt5 
     
To run the code use 
     python ./InfileGenerator.py
     
Developed by- Ajay Sagar
