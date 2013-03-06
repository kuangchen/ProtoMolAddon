#include <protomol/force/MetaForce.h>
#include <protomol/force/ForceGroup.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;
//____ MetaForce

void MetaForce::addToForceGroup(ForceGroup *forceGroup) {
  forceGroup->addMetaForce(this);
}

