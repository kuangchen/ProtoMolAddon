%module ScalarStructure
%{
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;
%}

%include <protomol/type/Real.h>
%include <protomol/type/ScalarStructure.h>

%extend ProtoMol::ScalarStructure {
   void setTable(int index, Real value) {
      self->operator[](ProtoMol::ScalarStructure::Index(index)) = value;
   }

   Real getTable(int index) {
      return self->operator[](ProtoMol::ScalarStructure::Index(index));
   }
}
