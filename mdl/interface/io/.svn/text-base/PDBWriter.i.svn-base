%module PDBWriter
%{
#include <protomol/io/PDBWriter.h>
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;
%}

%include <protomol/type/Real.h>
%include "std_string.i"
%include "std_vector.i"
%template() std::vector<PDB::Atom>;
%template() std::vector<PDB::Ter>;
%include <interface/type/Vector3DBlock.i>
%include <protomol/type/PDB.h>
%include <protomol/io/File.h>
%include <protomol/io/Writer.h>
%include <protomol/io/PDBWriter.h>
