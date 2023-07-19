#!/bin/bash
sed -i 's/gnome-terminal/dbus-launch gnome-terminal/g' MicroSim.py
sed -i 's/gnome-terminal/dbus-launch gnome-terminal/g' ./resources/functions/linux.py
sed -i 's/gnome-terminal/dbus-launch gnome-terminal/g' ./resources/functions/windows.py
sed -i 's/mpirun.mpich/mpirun/g' ./resources/functions/linux.py
sed -i 's/mpirun.mpich/mpirun/g' ./resources/functions/windows.py
mkdir bin

