%module XYZWriter
%{
#include <protomol/io/XYZWriter.h>
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;
%}

%include <protomol/type/Real.h>
%include "std_string.i"
%include "std_vector.i"
%template() std::vector<std::string>;
%template() std::vector<ProtoMol::Atom>;
%template() std::vector<ProtoMol::AtomType>;
%include <protomol/type/XYZ.h>
%include <interface/type/Vector3DBlock.i>
%include <protomol/topology/Atom.h>
%include <protomol/topology/AtomType.h>
%include <protomol/io/File.h>
%include <protomol/io/Writer.h>
%include <protomol/io/XYZWriter.h>
