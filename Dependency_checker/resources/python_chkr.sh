#!/bin/sh

#checking for python

python_version ()
{
python_version_required="$1"

echo "Checking for python compiler"
echo ""
python3 --version | grep Python >> version_log.txt
python_version_status=$?
if [ "$python_version_status" -eq 0 ]; then
   echo "exit-status is $python_version_status so Python exists \e[1;32m STATUS : OK \e[0m"
   python3 --version | grep Python
else
   echo "python 3 doesnt exist. Please install and try again \e[1;31m STATUS : NOT FOUND \e[0m "
fi  
#echo >> version_log.txt



python_ver_req=`echo $python_version_required | sed 's/\.//g'`
py_v_r=`echo "$python_ver_req/10" | bc`
python_ver_current=`python3 --version | grep -oP '(?<=Python )[^ ]*' | sed 's/\.//g'`
cond=`echo "$python_ver_current/10 >= $py_v_r" | bc`

if [ "$cond" -eq 1 ]; then
      echo "Python Version:`python3 --version | grep -oP '(?<=Python )[^ ]*' ` is >= $python_version_required  \e[0;32m STATUS : OK \e[0m  "
   else
      echo "Please install Python Version >= ($python_version_required) and try again. \e[1;36m STATUS : WRONG VERSION \e[0m "
     
    fi
}
