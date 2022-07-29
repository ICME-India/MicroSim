#!/bin/sh


pandas_chkr ()
{

echo "**************************************************************"
echo "Checking for pandas package"
echo ""
pip3 show pandas
pandas_curr_ver=`pip3 show pandas | grep Version`
pandas_status=$?
if [ "$pandas_status" -eq 0 ]; then
   echo ""
   echo "exit-status is $sympy_status so Pandas exists \e[1;32m STATUS : OK \e[0m "   
else   
   echo "Pandas doesnt exist please install the package and try again \e[1;31m STATUS : NOT FOUND \e[0m  "
   #exit
fi 
echo "pandas Current $pandas_curr_ver"
echo Pandas $pandas_curr_ver >> version_log.txt
echo >> version_log.txt
}
