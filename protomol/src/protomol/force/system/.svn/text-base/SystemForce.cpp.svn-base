#include <protomol/force/system/SystemForce.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/force/system/SystemCompareForce.h>
#include <protomol/force/system/SystemTimeForce.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;
//____ SystemForce

void SystemForce::addToForceGroup(ForceGroup *forceGroup) {
  forceGroup->addSystemForce(this);
}

CompareForce *SystemForce::makeCompareForce(Force *actualForce,
                                            CompareForce *compareForce) const {
  return new SystemCompareForce(actualForce, compareForce);
}

TimeForce *SystemForce::makeTimeForce(Force *actualForce) const {
  return new SystemTimeForce(actualForce);
}

