#!/bin/bash 

py_vt()

 {

echo "**************************************************************"
echo "Checking for vtk package"
echo ""
pip3 show vtk 
vtk_curr_ver=`pip3 show vtk | grep Version`
vtk_status=$?
if [ "$vtk_status" -eq 0 ]; then
   echo ""
   echo "exit-status is $vtk_status so vtk exists \e[1;32m STATUS : OK \e[0m "   
else   
   echo "vtk doesnt exist please install the package and try again \e[1;31m STATUS : NOT FOUND \e[0m  "
   #exit
fi 
echo "vtk Current $vtk_curr_ver"
echo vtk $vtk_curr_ver >> version_log.txt
echo >> version_log.txt
}
