#!/bin/sh


CL_Devices ()
{
echo "**************************************************************"
echo "checking for OpenCl installation"

cd OpenCL_DeviceInfo
sh compile.sh >/dev/null 2>&1
sh run.sh 
#sh run.sh | grep "Number of" &> /dev/null



initialize_devs ()
{
#echo "successul calling function"
no_of_devs=0
#no_of_platfs=0
}


initialize_platfs ()
{
#echo "successul calling function"
#no_of_devs=0
no_of_platfs=0
}


no_of_devs=`sh run.sh | grep -oP '(?<=Number of devices: )[^ ]*'`
stat=$?
#echo "Status of no of devices $stat"
if [ "$stat" -ge 1 ]; then 
 initialize_devs
 #no_of_devs=0
 #echo "No of Devices = $no_of_devs"
 #echo "Status of no of devices $stat"
fi 


 no_of_platfs=`sh run.sh | grep -oP '(?<=Number of OpenCL platforms: )[^ ]*'`
 stat1=$?
 #echo "Status of no of plats $stat1"
if [ "$stat1" -ge 1 ]; then
 initialize_platfs
 #echo "HIIIIIIIIIIIIIIIIIII  Status of no of plats $stat1"
 #echo "brooooooooooooooooooooo    No of platforms = $no_of_plats"
fi
#if [ "$no_of_devs" -eq 0 ]; then
if [ "$no_of_devs" -gt 0 -a "$no_of_platfs" -gt 0 ]; then
    echo " Number of devices = $no_of_devs. \n Number of Platforms= $no_of_platfs \n OpenCl is installed and properly configured. \e[1;32m STATUS : OK \e[0m"
else
    echo "Number of devices = $no_of_devs. \n Number of Platforms= $no_of_platfs \n OpenCl is not installed or it isnt properly configured  \n. \e[1;31m STATUS : NOT FOUND \e[0m"
fi
cd ..
}
