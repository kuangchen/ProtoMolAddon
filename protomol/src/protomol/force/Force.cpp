#include <protomol/force/Force.h>

#include <protomol/force/CompareForce.h>
#include <protomol/base/Exception.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ Force
const string Force::scope("Force");

Force *Force::make(const vector<Value> &values) const {
  assertParameters(values);
  return adjustAlias(doMake(values));
}

CompareForce *Force::makeCompareForce(Force *actualForce,
                                      CompareForce *compareForce) const {
  report << error << "No support to compare " <<
  (actualForce ? actualForce->getId() :
   string("null pointer")) << " vs. " <<
  (compareForce ? compareForce->getId() :
   string("null pointer")) << "." <<
  endr;
  return NULL;
}

TimeForce *Force::makeTimeForce(Force *actualForce) const {
  report << error << "No support to time " <<
  (actualForce ? actualForce->getId() :
   string("null pointer")) << "." <<
  endr;
  return NULL;
}

void Force::setParameters(vector<Value> values) {
  assertParameters(values);
  doSetParameters(values);
}

