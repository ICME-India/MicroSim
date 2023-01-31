### Requirements

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

To view the plots, gnuplot can be used, which can be installed using the following command in the terminal:

> sudo apt install gnuplot-x11

It can be launched using the command:

> gnuplot

It is advised to refer to the OpenFOAM documentions to understand the methods involved while using. The below link can be useful:

* https://cfd.direct/openfoam/documentation/


### OpenFOAM modules

Copy the modules to the OpenFOAM run directory. For instance, to copy the directory PFBinary:

> cp -r PFBinary $FOAM_RUN

Solver has to be compiled from the corresponding solver directory. For instance:

> cd $FOAM_RUN/PFBinary/PFBinary

> wclean

> wmake

Finally, switch to the case directory, e.g. coolingAlZn:

> cd $FOAM_RUN/PFBinary/coolingAlZn

To run the case with the default parameters, execute the following:

> ./Allclean

> ./Allrun

Note: Allclean and Allrun must be set as executables

To check the results in ParaView:

> paraview binary.foam

