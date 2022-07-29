#!/bin/sh

nvcc_checker()
{
 echo "**************************************************************"
echo "Checking for Nvidia CUDA Compiler package"
 echo "**************************************************************"
echo ""
echo "Nvidia CUDA Compiler's Required version: $1"   
 nvcc -v | grep -oP [?<=release]   
 nvcc_status=$?

#########################################################
 nvcc_ver_current=`nvcc -v | grep -oP '(?<=release )[^ ]*' | sed 's/\,//g'`
#echo $nvcc_status
#echo $nvcc_ver_current
 echo $nvcc_ver_req=$1
if [ "$nvcc_status" -eq 0 ]; then
cond3=`echo "$nvcc_ver_current >= $nvcc_ver_reqd" | bc`
    if [ "$cond3" -eq 1 ]; then
         echo "The required nvcc Version: $1 exists \e[0;32m STATUS : OK \e[0m"
          echo nvcc Version : $nvcc_version_current >> version_log.txt
    else
      echo "nvcc version : "$nvcc_version_current" exists but the required version is "$1" Please install and try again. \e[1;36m STATUS : WRONG VERSION \e[0m"
    fi
else
   echo "nvcc doesnt exist please install version "$1" and try again  \e[0;31m STATUS : NOT FOUND \e[0m  "
   
fi 

echo >> version_log.txt 
}
