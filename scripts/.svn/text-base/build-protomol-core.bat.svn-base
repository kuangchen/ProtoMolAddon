@echo off

REM assume check out is done

IF NOT %BUILD_ROOT%x == x GOTO CD_BUILD_ROOT
echo Please set BUILD_ROOT
GOTO END

:CD_BUILD_ROOT
cd %BUILD_ROOT%
IF %ERRORLEVEL% EQU 0 GOTO CHECK_BOOST
echo %BUILD_ROOT% does not exist
GOTO END

:CHECK_BOOST
IF EXIST %BOOST_HOME% GOTO BUILD_ROOT_SET
echo BOOST_HOME not set or invalid
GOTO END

:BUILD_ROOT_SET
set SCONSOPTS=debug=0 optimize=1 fah=1 simtk_lapack=1
set LIBBZ2_HOME=%BUILD_ROOT%\fah\libbzip2
set LIBFAH_HOME=%BUILD_ROOT%\fah\libfah
set PROTOMOL_HOME=%BUILD_ROOT%\protomol\protomol
set TARGET=

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
set DIR=%BUILD_ROOT%\protomol\fah-core
goto DO_PART
:PROTOMOL_CORE_DONE

set NAME=MSVC_PROTOMOL
set DIR=%PROTOMOL_HOME%
set TARGET=ProtoMol.vcproj
goto DO_PART
:MSVC_PROTOMOL_DONE

GOTO END

REM **************************************************************************
:DO_PART
IF EXIST %NAME%.done GOTO PART_DONE

echo Building %NAME%

cd %DIR%
IF %ERRORLEVEL% EQU 0 GOTO PART_BUILD
echo ERROR changing to directory %DIR%
GOTO END

:PART_BUILD
CALL scons %SCONSOPTS% %TARGET% %1
IF %ERRORLEVEL% EQU 0 GOTO PART_DONE
echo ERROR building %NAME%
GOTO END

:PART_DONE
cd %BUILD_ROOT%
echo done > %NAME%.done
echo %NAME% done
goto %NAME%_DONE


REM **************************************************************************
:END
cd %BUILD_ROOT%
