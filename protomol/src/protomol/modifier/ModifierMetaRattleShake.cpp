#include <protomol/modifier/ModifierMetaRattleShake.h>
#include <protomol/topology/Topology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/ProtoMolApp.h>

using namespace std;
using namespace ProtoMol;

//____ ModifierMetaRattleShake
ModifierMetaRattleShake::ModifierMetaRattleShake(Real eps, int maxIter,
                                                 bool all, int order) :
  Modifier(order), myEpsilon(eps), myMaxIter(maxIter), myAll(all), myListOfConstraints(0) {}

void ModifierMetaRattleShake::doInitialize() {
  myLastPositions = app->positions;

  // ... maybe it's a second inittialize, add the old back.
  app->topology->degreesOfFreedom +=
    app->topology->bondRattleShakeConstraints.size();

  buildRattleShakeBondConstraintList(app->topology,
                                     app->topology->bondRattleShakeConstraints, myAll);
  // this list contains bonded pairs, and UB-bonded pairs excluding
  // (heavy atom)-H pairs and (heavy)-(heavy) pairs

  // subtract the # of constraints from the # of degrees of freedom of the 
  // system. This is needed so that we get the correct temperature
  app->topology->degreesOfFreedom -=
    app->topology->bondRattleShakeConstraints.size();

  myListOfConstraints = &(app->topology->bondRattleShakeConstraints);
}

void ModifierMetaRattleShake::
getParameters(vector<Parameter> &parameters) const {
  parameters.push_back
    (Parameter("-epsilon",
               Value(myEpsilon, ConstraintValueType::NotNegative())));
  parameters.push_back
    (Parameter("-maxIter",
               Value(myMaxIter, ConstraintValueType::NotNegative())));
  parameters.push_back
    (Parameter("-all",
               Value(myAll, ConstraintValueType::NoConstraints())));
}
