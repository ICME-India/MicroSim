#!/bin/sh


h5pcc_chkr ()
{
echo "Checking for h5pcc compiler"
echo ""
h5pcc --version
h5pcc_status=$?
 #shouldnot use echo as it takes echo as a success
if [ "$h5pcc_status" -eq 0 ]; then
   echo h5pcc Version:`h5pcc --version | grep gcc | grep .` >> version_log.txt
   echo "exit-status is $h5pcc_status so h5pcc exists  \e[1;32m STATUS : OK \e[0m"
   
else
   echo "h5pcc isnt installed or configured properly, please install and try again \e[0;31m STATUS : NOT FOUND \e[0m "
   
fi   
echo >> version_log.txt
}
