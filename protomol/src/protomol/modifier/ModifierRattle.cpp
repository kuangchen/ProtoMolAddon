#include <protomol/modifier/ModifierRattle.h>
#include <protomol/integrator/Integrator.h>
#include <protomol/topology/Topology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/PMConstants.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ ModifierRattle
ModifierRattle::ModifierRattle() : ModifierMetaRattle(0, 0, true, 0) {}
ModifierRattle::ModifierRattle(Real eps, int maxIter, bool all, int order) :
  ModifierMetaRattle(eps, maxIter, all, order) {}

void ModifierRattle::doExecute(Integrator *i) {
  // estimate the current error in all velocity constraints
  Real error = calcError();

  // delta_t
  Real dt = i->getTimestep() / Constant::TIMEFACTOR;

  int iter = 0;
  while (error > myEpsilon) {
    for (unsigned int i = 0; i < myListOfConstraints->size(); i++) {
      // find the ID#s of the two atoms in the current constraint
      int a1 = (*myListOfConstraints)[i].atom1;
      int a2 = (*myListOfConstraints)[i].atom2;

      // reciprocal atomic masses
      Real rM1 = 1 / app->topology->atoms[a1].scaledMass;
      Real rM2 = 1 / app->topology->atoms[a2].scaledMass;

      // now lets compute the lambdas.
      // compute the current bond vector
      Vector3D rab = app->positions[a1] - app->positions[a2];
      Real rabsq = rab.normSquared();

      // compute the current velocity vector
      Vector3D vab = app->velocities[a1] - app->velocities[a2];

      // dot product of distance and velocity vectors
      Real rvab = rab * vab;

      // compute the change in lambda
      Real gab = -rvab / (dt * (rM1 + rM2) * rabsq);
      Vector3D dp = rab * gab;

      // move the velocities based upon the multiplier
      app->velocities[a1] += dp * dt * rM1;
      app->velocities[a2] -= dp * dt * rM2;

      // the constraint adds a force to each atom since their positions
      // had to be changed.  This constraint force therefore contributes
      // to the atomic virial.  Note that the molecular virial is independent of
      // any intramolecular constraint forces.
      if (app->energies.virial()) app->energies.addVirial(dp * 2, rab);
    }

    // compute the error in all the velocity constraints after this RATTLE
    // iteration
    error = calcError();
    iter++;
    if (iter > myMaxIter) {
      report << warning << "maxIter = " << myMaxIter
             << " reached, but still not converged ... error is " << error <<
      endr;
      break;
    }
  }
}
