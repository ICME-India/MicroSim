#!/bin/sh

. ./resources/python_chkr.sh
. ./resources/pycalphad_pkgs.sh
. ./resources/sympy_pkg_chkr.sh
. ./resources/pandas_pkg_chkr.sh
. ./CL_Deviceinfo.sh
. ./resources/h5pcc_chkr.sh
. ./c_resources/c.sh


#	. ./C_Header_files/c_gsl_chkr.sh
h5pcc_chkr
python_version 3.7.0 #enter the required version number in place of 3.7.0
pycalphad_version 0.9.2 #enter the required version number in place of 0.9.2
sympy_chkr
pandas_chkr
gslchkr
fftw3chkr
CL_Devices

