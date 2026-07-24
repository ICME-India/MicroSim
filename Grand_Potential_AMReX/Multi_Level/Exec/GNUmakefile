# AMREX_HOME defines the directory in which we will find all the AMReX code.
AMREX_HOME ?= /opt/AMREX/amrex
#LIBRARY_LOCATIONS += /home/ext/apps/spack/opt/spack/linux-centos7-cascadelake/gcc-11.2.0/gsl-2.7-h5q52xdndmnqwjwnysbdt5c6ccbk4fv6/lib 
LIBRARIES += -lgsl -lgslcblas -lm

INCLUDE_LOCATIONS += ../Source #/home/ext/apps/spack/opt/spack/linux-centos7-cascadelake/gcc-11.2.0/gsl-2.7-h5q52xdndmnqwjwnysbdt5c6ccbk4fv6/include

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE
DEBUG      = FALSE

DIM        = 2
#DIM       = 3

COMP	   = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE 
USE_CUDA   = FALSE


Blocs   := . 

include Make.Adv
