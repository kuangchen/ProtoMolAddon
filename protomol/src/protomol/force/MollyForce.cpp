#include <protomol/force/MollyForce.h>
#include <protomol/force/ForceGroup.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;
//____ MollyForce

void MollyForce::addToForceGroup(ForceGroup *forceGroup) {
  forceGroup->addMollyForce(this);
}

