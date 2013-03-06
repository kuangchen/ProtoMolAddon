#include <protomol/modifier/ModifierMetaRattle.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/topology/Topology.h>

using namespace std;
using namespace ProtoMol;

//____ ModifierMetaRattle
ModifierMetaRattle::ModifierMetaRattle(Real eps, int maxIter, bool all, int order) :
  ModifierMetaRattleShake(eps, maxIter, all, order) {}

Real ModifierMetaRattle::calcError() const {
  // the error for the RATTLE algorithm is defined as
  // fabs( [v1-v2] * [r1-r2] ) < tolerance
  // which is the constraint imposed upon the velocities by RATTLE.
  // It is this constraint upon the velocities that allows us to compute
  // the multipliers (lambdas) at time t + delta_t

  Real error = 0;
  for (unsigned int i = 0; i < myListOfConstraints->size(); i++) {
    int a1 = (*myListOfConstraints)[i].atom1;
    int a2 = (*myListOfConstraints)[i].atom2;
    Vector3D vab = app->velocities[a1] - app->velocities[a2];
    Vector3D rab = app->positions[a1] - app->positions[a2];
    Real err = fabs(rab * vab);
    error += err;
  }

  return error /= myListOfConstraints->size();
}

