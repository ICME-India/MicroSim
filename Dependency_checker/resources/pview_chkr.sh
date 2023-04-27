#! /bin/sh


paraview_chkr()
{
paraview --version #| grep -oP '(?<=version) [^ ]*'
paraview_status=$?
if [ "$paraview_status" -eq 0 ]; then
    echo "exit-status is $paraview_status so Paraview exists \e[1;32m STATUS : OK \e[0m"
fi
}
