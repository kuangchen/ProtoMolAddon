#include "ModifierNPTShake.h"
#include <protomol/integrator/Integrator.h>
#include <protomol/topology/Topology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>

using namespace ProtoMol::Report;

namespace ProtoMol {


  //__________________________________________________ ModifierNPTShakeDetails
  ModifierNPTShakeDetails::ModifierNPTShakeDetails(Real eps, int maxIter, int order):ModifierMetaShake(eps,maxIter,true,order){}

  void ModifierNPTShakeDetails::doExecute(Integrator* myTheIntegrator){
    // estimate the current error in all bond constraints
    Real error = calcError();

    // delta_t
    Real dt = myTheIntegrator->getTimestep() / Constant::TIMEFACTOR;
    Real dtsq = dt * dt;

    // multiplicative constants for the box volume velocity and thermostat velocity
    const Real pVol_term = 0.5 * myTheIntegrator->getTimestep() * getEpsilonVel();
    const Real pEta_term = 0.5 * myTheIntegrator->getTimestep() * getEtaVel();
    const Real OneOverN_term = 1.0 + 1.0 / app->topology->atoms.size();

    int iter = 0;
    while(error > myEpsilon) {
      
      for(unsigned int i=0;i<myListOfConstraints->size();i++) {

        // find the ID#s of the two atoms in the current constraint
        int a1 = (*myListOfConstraints)[i].atom1;
        int a2 = (*myListOfConstraints)[i].atom2;

        // get the ID# of the molecule that this bond belongs to
        int Mol = app->topology->atoms[a1].molecule;

        // reciprocal atomic masses
        Real rM1 = 1/app->topology->atoms[a1].scaledMass;
        Real rM2 = 1/app->topology->atoms[a2].scaledMass;
        Real m1_over_M = 1 / (rM1 * app->topology->molecules[Mol].mass);
        Real m2_over_M = 1 / (rM2 * app->topology->molecules[Mol].mass);

        // multiplicative constants due to change in box volume
        // (for constant pressure simulations -- see Equation 56 of G. Kalibaeva,
        //  M. Ferrario, and G. Ciccotti, "Constant pressure-constant temperature molecular
        //  dynamics: a correct constrained NPT ensemble using the molecular virial", Mol. Phys.
        //  101(6), 2003, p. 765-778)
        Real Exp1_atom1 = exp(pVol_term * m1_over_M);
        Real Exp1_atom2 = exp(pVol_term * m2_over_M);
        Real Exp2_atom1 = (1 + pVol_term*OneOverN_term*m1_over_M) * exp( -pEta_term - pVol_term*OneOverN_term*m1_over_M);
        Real Exp2_atom2 = (1 + pVol_term*OneOverN_term*m2_over_M) * exp( -pEta_term - pVol_term*OneOverN_term*m2_over_M);

        // get the target bond distance for this constraint
        Real restLength = (*myListOfConstraints)[i].restLength;
        
        // now lets compute the lambdas.
        // compute the current bond vector
        Vector3D pab = (app->positions)[a1] - (app->positions)[a2];   
        Real pabsq = pab.normSquared();
        Real rabsq = restLength*restLength;
        
        // compute the difference between the target bond length and
        // the actual bond length
        Real diffsq = rabsq - pabsq; //-g^k()

        // compute the bond vector from the previous timestep
        Vector3D rab = myLastPositions[a1] - myLastPositions[a2];
        Real rpab = rab * pab;

        // calculate the constraint force, or multiplier
        Real gab = diffsq / (2 * dtsq * (rM1 * Exp1_atom1 + rM2 * Exp1_atom2) * rpab);
        Vector3D dp = rab * gab;

        // move the positions based upon the multiplier
        (app->positions)[a1] += dp * dtsq * rM1 * Exp1_atom1;
        (app->positions)[a2] -= dp * dtsq * rM2 * Exp1_atom2;

        // move the velocities based upon the multiplier
        (app->velocities)[a1] += dp * dt * rM1 * Exp2_atom1;
        (app->velocities)[a2] -= dp * dt * rM2 * Exp2_atom2;

        // the constraint adds a force to each atom since their positions
        // had to be changed.  This constraint force therefore contributes
        // to the atomic virial.  Note that the molecular virial is independent of
        // any intramolecular constraint forces.
        if (app->energies.virial()) app->energies.addVirial(dp*2,rab);
      }
     
      // compute the error in all the bond constraints after this SHAKE iteration
      error = calcError();
      iter ++;
      if(iter > myMaxIter) {
	report << warning << "Shake maxIter = " << myMaxIter 
	       << " reached, but still not converged ... error is "<<error<<endr;
	break;
      }      
    }

    // store the old positions
    myLastPositions = (app->positions);
  }

}
