#include <protomol/force/extended/ExtendedForce.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/force/extended/ExtendedCompareForce.h>
#include <protomol/force/extended/ExtendedTimeForce.h>

using namespace ProtoMol;
//____ ExtendedForce

void ExtendedForce::addToForceGroup(ForceGroup *forceGroup) {
  forceGroup->addExtendedForce(this);
}

CompareForce *ExtendedForce::makeCompareForce(Force *actualForce,
                                              CompareForce *compareForce)
const {
  return new ExtendedCompareForce(actualForce, compareForce);
}

TimeForce *ExtendedForce::makeTimeForce(Force *actualForce) const {
  return new ExtendedTimeForce(actualForce);
}

