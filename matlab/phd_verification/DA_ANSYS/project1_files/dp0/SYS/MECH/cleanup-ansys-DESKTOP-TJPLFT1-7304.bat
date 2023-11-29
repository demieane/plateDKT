@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="DESKTOP-TJPLFT1" (taskkill /f /pid 15532)
if /i "%LOCALHOST%"=="DESKTOP-TJPLFT1" (taskkill /f /pid 7304)

del /F cleanup-ansys-DESKTOP-TJPLFT1-7304.bat
