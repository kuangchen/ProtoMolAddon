#include <protomol/type/GromacsParameters.h>

using namespace ProtoMol;

void GromacsParameters::clear() {

   defaults.clear();
   bondTypes.clear();
   angleTypes.clear();
   dihedralTypes.clear();
   atomTypes.clear();
   gb_parameters.clear();
}
