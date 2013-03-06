@echo off

REM assume check out is done

IF NOT %BUILD_ROOT%x == x GOTO CD_BUILD_ROOT
echo Please set BUILD_ROOT
GOTO END

:CD_BUILD_ROOT
cd %BUILD_ROOT%
IF %ERRORLEVEL% EQU 0 GOTO SET_VARS
echo %BUILD_ROOT% does not exist
GOTO END

:SET_VARS
set LIBBZ2_HOME=%BUILD_ROOT%\fah\libbzip2
set LIBFAH_HOME=%BUILD_ROOT%\fah\libfah
set PROTOMOL_HOME=%BUILD_ROOT%\protomol\protomol

REM **************************************************************************
set NAME=LIBBZIP2
set DIR=%LIBBZ2_HOME%
goto DO_PART
:LIBBZIP2_DONE

set NAME=LIBFAH
set DIR=%LIBFAH_HOME%
goto DO_PART
:LIBFAH_DONE

set NAME=LIBPROTOMOL_FAH
set DIR=%PROTOMOL_HOME%
goto DO_PART
:LIBPROTOMOL_FAH_DONE

set NAME=PROTOMOL_CORE
set DIR=%BUILD_ROOT%\cores\protomol
goto DO_PART
:PROTOMOL_CORE_DONE

GOTO END


REM **************************************************************************
:DO_PART
IF EXIST %NAME%.done del %NAME%.done

echo Cleaning %NAME%

cd %DIR%
IF %ERRORLEVEL% EQU 0 GOTO PART_CLEAN
echo ERROR changing to directory %DIR%
GOTO END

:PART_CLEAN
CALL scons -c
IF %ERRORLEVEL% EQU 0 GOTO PART_DONE
echo ERROR cleaning %NAME%
GOTO END

:PART_DONE
cd %BUILD_ROOT%
GOTO %NAME%_DONE


REM **************************************************************************
:END
cd %BUILD_ROOT%
