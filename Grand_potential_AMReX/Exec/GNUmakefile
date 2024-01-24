# AMREX_HOME defines the directory in which we will find all the AMReX code.
AMREX_HOME ?= /opt/AMREX/amrex
#LIBRARY_LOCATIONS += /home/ext/apps/spack/opt/spack/linux-centos7-cascadelake/gcc-11.2.0/gsl-2.7-h5q52xdndmnqwjwnysbdt5c6ccbk4fv6/lib 


LIBRARIES += -lgsl -lgslcblas -lm

DEBUG        = FALSE
USE_MPI      = TRUE
USE_CUDA     = FALSE
PROFILE 	 = FALSE
USE_OMP      = FALSE
COMP         = gnu
DIM          = 2

#phasecount = 2
#compcount = 2
#DEFINES = -Dphasecount=$(phasecount)
#DEFINES = -Dcompcount=$(compcount)

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ../Source/Make.package
VPATH_LOCATIONS  += ../Source
INCLUDE_LOCATIONS += ../Source #/home/ext/apps/spack/opt/spack/linux-centos7-cascadelake/gcc-11.2.0/gsl-2.7-h5q52xdndmnqwjwnysbdt5c6ccbk4fv6/include


include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
