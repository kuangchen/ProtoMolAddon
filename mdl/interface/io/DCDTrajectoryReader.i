%module DCDTrajectoryReader
%{
#include <protomol/io/DCDTrajectoryReader.h>
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;
%}

%include <protomol/type/Real.h>
%include <protomol/type/XYZ.h>
%include "std_string.i"
%include "std_vector.i"
%include <protomol/io/DCDTrajectoryReader.h>


%extend ProtoMol::DCDTrajectoryReader {
   int numTrajectories() {return self->xyz->size();}
   int size() {
      return self->xyz->size();
   }
   Real getElement(int index2, int index3) {
      return (*(self->xyz))[index2][index3];
   }
}
