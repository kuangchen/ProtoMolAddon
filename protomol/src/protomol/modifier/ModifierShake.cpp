#include <protomol/modifier/ModifierShake.h>
#include <protomol/integrator/Integrator.h>
#include <protomol/topology/Topology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/PMConstants.h>

using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ ModifierShake
ModifierShake::ModifierShake() : ModifierMetaShake(0, 0, true, 0) {}
ModifierShake::ModifierShake(Real eps, int maxIter, bool all, int order) :
  ModifierMetaShake(eps, maxIter, all, order) {}

void ModifierShake::doExecute(Integrator *i) {
  // estimate the current error in all bond constraints
  Real error = calcError();

  // delta_t
  Real dt = i->getTimestep() / Constant::TIMEFACTOR;

  int iter = 0;
  while (error > myEpsilon) {
    for (unsigned int i = 0; i < myListOfConstraints->size(); i++) {
      // find the ID#s of the two atoms in the current constraint
      int a1 = (*myListOfConstraints)[i].atom1;
      int a2 = (*myListOfConstraints)[i].atom2;

      // get the target bond distance for this constraint
      Real restLength = (*myListOfConstraints)[i].restLength;

      // now lets compute the lambdas.
      Vector3D pab = app->positions[a1] - app->positions[a2];

      // compute the current bond vector
      Real pabsq = pab.normSquared();
      Real rabsq = restLength * restLength;

      // compute the difference between the target bond length and
      // the actual bond length
      Real diffsq = rabsq - pabsq;   //-g^k()

      // compute the bond vector from the previous timestep
      Vector3D rab = myLastPositions[a1] - myLastPositions[a2];
      Real rpab = rab * pab;

      // reciprocal atomic masses
      Real rM1 = 1 / app->topology->atoms[a1].scaledMass;
      Real rM2 = 1 / app->topology->atoms[a2].scaledMass;

      // calculate the constraint force, or multiplier
      Real gab = diffsq / (2 * (rM1 + rM2) * rpab);
      Vector3D dp = rab * gab;

      // move the positions based upon the multiplier
      app->positions[a1] += dp * rM1;
      app->positions[a2] -= dp * rM2;

      dp /= dt;

      // move the velocities based upon the multiplier
      app->velocities[a1] += dp * rM1;
      app->velocities[a2] -= dp * rM2;

      // the constraint adds a force to each atom since their positions
      // had to be changed.  This constraint force therefore contributes
      // to the atomic virial.  Note that the molecular virial is independent of
      // any intramolecular constraint forces.
      if (app->energies.virial()) {
        dp /= dt;
        app->energies.addVirial(dp * 2, rab);
      }
    }

    // compute the error in all the bond constraints after this SHAKE iteration
    error = calcError();
    iter++;
    if (iter > myMaxIter) {
      report << warning << "maxIter = " << myMaxIter
             << " reached, but still not converged ... error is "
             << error << endr;
      break;
    }
  }

  // store the old positions
  myLastPositions = app->positions;
}

