# AMREX_HOME defines the directory in which we will find all the AMReX code.
AMREX_HOME ?= /opt/AMREX/amrex
LIBRARY_LOCATIONS += /usr/lib
XTRALIBS += -lgsl -lgslcblas -lm

DEBUG        = FALSE
USE_MPI      = TRUE
USE_CUDA     = FALSE
USE_OMP      = FALSE
COMP         = gnu
DIM          = 2

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ../Source/Make.package
VPATH_LOCATIONS  += ../Source
INCLUDE_LOCATIONS += ../Source /usr/include

include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
