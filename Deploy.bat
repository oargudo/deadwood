:: x64 Release
cls

:: GL
XCOPY /Y "%GLEW_DIR%\bin\Release\x64\glew32.dll" "bin\x64\Debug\" 
XCOPY /Y "%GLEW_DIR%\bin\Release\x64\glew32.dll" "bin\x64\Release\" 

:: Deploy
%QT_DIR%\bin\windeployqt.exe --release bin\x64\Release\AppDeadwood.exe
%QT_DIR%\bin\windeployqt.exe bin\x64\Debug\AppDeadwood.exe