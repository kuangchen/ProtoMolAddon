#include <protomol/modifier/ModifierMetaShake.h>
#include <protomol/topology/Topology.h>
#include <protomol/ProtoMolApp.h>

using namespace ProtoMol;

//____ ModifierMetaShake
ModifierMetaShake::ModifierMetaShake(Real eps, int maxIter, bool all, int order) :
  ModifierMetaRattleShake(eps, maxIter, all, order) {}

Real ModifierMetaShake::calcError() const {
  // the error is defined as < fabs(dist - restLength)/restLength >,
  // which is approximated as 
  // < fabs(dist^2 - restLength^2)/(2*restLength*restLength >.
  // This approximation avoids the sqrt operation to compute the distance.
  // cf. Krautler, van Gunsteren, et al J. Comput. Chem., 22(5) 501--508 (2001)
  // "a fast shake algorithm to solve distance constraint equations for small
  // molecules in molecular dynamics simulations"
  Real error = 0;
  for (unsigned int i = 0; i < myListOfConstraints->size(); i++) {
    int a1 = (*myListOfConstraints)[i].atom1;
    int a2 = (*myListOfConstraints)[i].atom2;
    Real restLengthSquared = power<2>((*myListOfConstraints)[i].restLength);
    Real err = fabs(
      (app->positions[a2] -
       app->positions[a1]).normSquared() -
      restLengthSquared) / (2.0 * restLengthSquared);
    error += err;
  }

  return error /= myListOfConstraints->size();
}
