#!/bin/sh

sympy_chkr ()
{
echo "**************************************************************"
echo "Checking for sympy package"
echo ""
pip3 show sympy
Sympy_curr_ver=`pip3 show sympy | grep Version`
sympy_status=$?
if [ "$sympy_status" -eq 0 ]; then
   echo "exit-status is $sympy_status so Sympy exists                \e[0;32m STATUS : OK \e[0m "   
else
   echo "Sympy doesnt exist please install the package and try again \e[0;31m STATUS : NOT FOUND \e[0m "
   #exit 
fi 
echo "Sympy Current $Sympy_curr_ver"
echo Sympy $Sympy_curr_ver >> version_log.txt
echo >> version_log.txt
}
