#!/bin/bash
{
mkdir /opt/MicroSim
chmod +x bin/MicroSim
cp -r bin/. /opt/MicroSim/
cp /opt/MicroSim/MicroSim.desktop /usr/share/applications/MicroSim.desktop
echo "source /opt/MicroSim/MSbashrc" >> ~/.bashrc
mkdir -p /home/$SUDO_USER/MicroSim/bin
cp -r solvers/. /home/$SUDO_USER/MicroSim/
chown -R $SUDO_USER:$SUDO_USER /home/$SUDO_USER/MicroSim
} || {
echo "Oops ! Installation failed due to some error. Always run this script as a super user (i.e., sudo sh install.sh)"
exit
}
echo "--------------------------------------------------------------------"
echo ""
echo "             __  __ _               ____  _           
            |  \/  (_) ___ _ __ ___/ ___|(_)_ __  ___  
            | |\/| | |/ __| '__/ _ \___ \| | '_  \` _ \ 
            | |  | | | (__| | | (_) |__) | | | | | | |
            |_|  |_|_|\___|_|  \___/____/|_|_| |_| |_|
            
          M I C R O S T R U C T U R E   S I M U L A T O R
"
echo ""
echo "--------------------------------------------------------------------"
echo ""
echo ""
echo " MicroSim installed successfully..."
echo ""
echo ""

