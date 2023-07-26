#!/bin/bash

sed -i 's/dbus-launch gnome-terminal/gnome-terminal/g' MicroSim.py
sed -i 's/dbus-launch gnome-terminal/gnome-terminal/g' ./resources/functions/linux.py
sed -i 's/dbus-launch gnome-terminal/gnome-terminal/g' ./resources/functions/windows.py
sed -i 's/mpirun/mpirun.mpich/g' ./resources/functions/linux.py
sed -i 's/mpirun.mpich/mpirun/g' ./resources/functions/windows.py
rm -r bin

