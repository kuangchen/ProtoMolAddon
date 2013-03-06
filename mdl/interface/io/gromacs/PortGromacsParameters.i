%module PortGromacsParameters
%{
#include <protomol/io/gromacs/PortGromacsParameters.h>
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;
%}

%include <std_string.i>
//%include <protomol/type/GromacsParameters.h>
//%include <protomol/type/GromacsTopology.h>
//%include <protomol/type/PAR.h>
//%include <protomol/type/PSF.h>
%include <protomol/io/gromacs/PortGromacsParameters.h>
