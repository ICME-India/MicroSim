#! /bin/sh


gnuplot_chkr()
{
    echo "***********************************************************"
    echo"Checking for gnuplot"
    echo "***********************************************************"
gnuplot --version #| grep -oP '(?<=gnuplot ) [^ ]*'
gnuplot_status=$?
if [ "$gnuplot_status" -eq 0 ]; then
    echo "exit-status is $gnuplot_status so gnuplot exists \e[1;32m STATUS : OK \e[0m"
else 
    echo  "\e[0;31m STATUS : NOT FOUND \e[0m "    
fi
}
