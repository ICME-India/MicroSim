# AMREX_HOME defines the directory in which we will find all the AMReX code.
AMREX_HOME ?= /opt/AMREX/amrex
#LIBRARY_LOCATIONS += /opt/homebrew/lib/
LIBRARIES += -lgsl -lgslcblas -lm

INCLUDE_LOCATIONS += ../Source #/opt/homebrew/include/

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
