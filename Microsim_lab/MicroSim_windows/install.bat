mkdir "C:\Program Files\MicroSim"
Xcopy /E /H /C /I "%~dp0\MicroSim\" "C:\Program Files\MicroSim\"
Xcopy "C:\Program Files\MicroSim\MicroSim.lnk" "%USERPROFILE%\Desktop\"
Xcopy "C:\Program Files\MicroSim\MicroSim.lnk" "%AppData%\Microsoft\Windows\Start Menu\Programs\"
Xcopy "%~dp0\solvers\" "%userprofile%\MicroSim\Solvers\" /d /e
pause
