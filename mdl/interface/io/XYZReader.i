%module XYZReader
%{
#include <protomol/io/XYZReader.h>
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;
%}

%include <protomol/type/Real.h>
%include "std_string.i"
%include <protomol/type/Vector3DBlock.h>
%include <protomol/type/XYZ.h>
%include <protomol/io/File.h>
%include <protomol/io/Reader.h>
%include <protomol/io/XYZReader.h>
