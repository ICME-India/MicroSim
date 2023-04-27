#!/bin/sh


#Searching for python Packages

pycalphad_version ()
{

echo "Checking for pycalphad package"
echo ""
echo "Required version: $1"
pycalphad_ver_req=`echo $1 | sed 's/\.//g'`
pycalphad_version_current=`pip3 show pycalphad | grep -oP '(?<=Version: )[^ ]*'` 
pycalphad_status=$?
 pycalphad_ver_current=`pip3 show pycalphad | grep -oP '(?<=Version: )[^ ]*' | sed 's/\.//g'`
#echo $pycalphad_status
#echo $pycalphad_ver_current
 #echo $pycalphad_ver_req
if [ "$pycalphad_status" -eq 0 ]; then
cond2=`echo "$pycalphad_ver_current == $pycalphad_ver_req" | bc`
 if [ "$cond2" -eq 1 ]; then
      echo "The required Pycalphad Version: $PYCALPHAD_VERSION_REQUIRED exists \e[0;32m STATUS : OK \e[0m"
      echo Pycalphad Version : $pycalphad_version_current >> version_log.txt
 else
      echo "Pycalphad version : "$pycalphad_version_current" exists but the required version is "$PYCALPHAD_VERSION_REQUIRED" Please install and try again. \e[1;36m STATUS : WRONG VERSION \e[0m"
  fi
else
   echo "Pycalphad doesnt exist please install version 0.9.2 and try again  \e[0;31m STATUS : NOT FOUND \e[0m  "
   
fi 

echo >> version_log.txt


#echo "Current version: $pycalphad_version_current"
}
