@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="DESKTOP-TJPLFT1" (taskkill /f /pid 16724)
if /i "%LOCALHOST%"=="DESKTOP-TJPLFT1" (taskkill /f /pid 20032)

del /F cleanup-ansys-DESKTOP-TJPLFT1-20032.bat
