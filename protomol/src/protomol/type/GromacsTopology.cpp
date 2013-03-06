#include <protomol/type/GromacsTopology.h>

using namespace ProtoMol;

void GromacsTopology::clear() {

   atoms.clear();
   bonds.clear();
   pairs.clear();
   angles.clear();
   rb_dihedrals.clear();
   dihedrals.clear();
}
