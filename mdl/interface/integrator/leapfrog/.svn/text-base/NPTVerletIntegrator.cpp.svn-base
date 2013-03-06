//  -----------------------------------------------------------------------  //
//  explicit, time-reversible integrator for NPT dynamics                    //
//                                                                           //
//  Unless modified, this integrator uses the molecular virial to control    //
//  the pressure. -- TIM                                                     //
//  -----------------------------------------------------------------------  //

#include "NPTVerletIntegrator.h"
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/base/PMConstants.h>
#include <protomol/topology/TopologyUtilities.h>
#include "ModifierPreForceThermostat.h"
#include "ModifierPreForceBarostat.h"
#include "ModifierPostForceThermostat.h"
#include "ModifierPostForceBarostat.h"
#include "ModifierNPTShake.h"
#include "ModifierNPTRattle.h"
#include <protomol/ProtoMolApp.h>

using std::vector;
using std::string;

namespace ProtoMol {

  //  -----------------------------------------------------------------------  //
  //  Keyword.
  const string NPTVerletIntegrator::keyword("NPTVerlet");

  //  -----------------------------------------------------------------------  //
  //  Default or empty constructor
  NPTVerletIntegrator::NPTVerletIntegrator(): STSIntegrator(),
                                              myTargetTemp(0.0),
                                              myTargetPres(0.0),
                                              myTauT(0.0),
                                              myTauV(0.0),
					      myTauP(0.0),
                                              kbT(0.0) {}

  //  -----------------------------------------------------------------------  //
  //  Constructor
  NPTVerletIntegrator::NPTVerletIntegrator(Real timestep,
					   Real temperature,
					   Real pressure,
					   Real tauT,
					   Real tauV,
					   Real tauP,
					   ForceGroup *overloadedForces)
    : STSIntegrator(timestep, overloadedForces),
      myTargetTemp(temperature), 
      myTargetPres(pressure),
      myTauT(tauT),
      myTauV(tauV),
      myTauP(tauP),
      kbT(temperature * Constant::BOLTZMANN){}

  //  -----------------------------------------------------------------------  //
  //  Thermostat -- prior to force calculations
  void NPTVerletIntegrator::PreForceThermostat() {

    //  Timestep.  Units: (fs)
    const Real halfDeltaT = 0.5 * getTimestep();   
    //  Twice the kinetic energy.  Units: (kcal / mol)
    const Real twiceKE = 2. * kineticEnergy(app->topology, &app->velocities);


    //  Advance the particle thermostat variable.
    myEta += myEtaVel * halfDeltaT;
    //  Advance the particle thermostat variable velocity.
    myEtaVel += (twiceKE - app->topology->degreesOfFreedom * kbT) * halfDeltaT / Qo;
  }

  //  -----------------------------------------------------------------------  //
  //  Thermostat -- after force calculations
  void NPTVerletIntegrator::PostForceThermostat() {

    //  Timestep.  Units: (fs)
    const Real halfDeltaT = 0.5 * getTimestep();    
    //  Get the new KE.
    const Real twiceKE = 2. * kineticEnergy(app->topology, &app->velocities);

       
    //  Advance the particle thermostat variable velocity.
    myEtaVel += (twiceKE - app->topology->degreesOfFreedom * kbT) * halfDeltaT / Qo;
    //  Advance the particle thermostat variable.
    myEta += myEtaVel * halfDeltaT;
   
    //  New particle thermostat kinetic energy.
    Real pEta = (Qo * 0.5) * (myEtaVel * myEtaVel);
    //  New particle thermostat potential energy.
    Real VEta = myEta * app->topology->degreesOfFreedom * kbT;

    //  Add the energy from the extended system thermostat.
    (app->energies)[ScalarStructure::INTEGRATOR] += pEta + VEta;
  }

  //  -----------------------------------------------------------------------  //
  //  Barostat -- prior to force calculations
  void NPTVerletIntegrator::PreForceBarostat() {

    //  Timestep.  Units: (fs)
    const Real halfDeltaT = 0.5 * getTimestep();    
    //  Twice the kinetic energy.  Units: (kcal / mol)
    const Real twiceKE = 2. * kineticEnergy(app->topology, &app->velocities);
    const Real MolKE = molecularKineticEnergy(app->topology, &app->velocities);


    //  Advance the volume thermostat variable.
    myEtaV += myEtaVolVel * halfDeltaT;
    //  Advance the volume thermostat variable velocity.
    myEtaVolVel += (W * myEpsilonVel * myEpsilonVel - kbT) * halfDeltaT / Qv;
    
    //  Calculate the current molecular pressure.  Units: (bar)
    Real currentPres = computeMolecularPressure(&app->energies, myVolume, MolKE);

    //  Advance the box volume velocity.  Divide by PRESSUREFACTOR to get
    //  correct (fs)^-1 units.
    myEpsilonVel += (3.0 * (myVolume) * (currentPres - myTargetPres)) *
      halfDeltaT / (Constant::PRESSUREFACTOR * W);
    //  Advance the box volume velocity due to thermostat forces.
    myEpsilonVel += (3.0 * twiceKE / app->topology->degreesOfFreedom) * halfDeltaT / W;
    myEpsilonVel *= exp(-myEtaVolVel * halfDeltaT);

    //  Advance the box volume.
    Real fac = exp(3.0 * myEpsilonVel * 2. * halfDeltaT);
    myVolume *= fac;
    app->topology->rescaleVolume(fac);
    uncache();
  }

  //  -----------------------------------------------------------------------  //
  //  Barostat -- after force calculations
  void NPTVerletIntegrator::PostForceBarostat() {

    //  Timestep.  Units: (fs)
    const Real halfDeltaT = 0.5 * getTimestep();   
    //  Get the new KE.
    const Real twiceKE = 2. * kineticEnergy(app->topology, &app->velocities);
    const Real MolKE = molecularKineticEnergy(app->topology, &app->velocities);

    //  Calculate the current pressure.  Units: (bar)
    Real currentPres = computeMolecularPressure(&app->energies, myVolume, MolKE);

    //  Advance the box volume velocity.
    myEpsilonVel *= exp(-myEtaVolVel * halfDeltaT);
    myEpsilonVel += (3.0 * twiceKE / app->topology->degreesOfFreedom) * halfDeltaT / W;
    //  Advance the box volume velocity.  Divide by PRESSUREFACTOR to get
    //  correct (fs)^-1 units.
    myEpsilonVel += (3.0 * (myVolume) * (currentPres - myTargetPres)) *
      halfDeltaT / (Constant::PRESSUREFACTOR * W);

    //  Advance the volume thermostat variable velocity.
    myEtaVolVel += (W * myEpsilonVel * myEpsilonVel - kbT) * halfDeltaT / Qv;
    //  Advance the volume thermostat variable.
    myEtaV += myEtaVolVel * halfDeltaT;


    //  New volume thermostat kinetic energy.
    Real pEtaV = (Qv * 0.5) * (myEtaVolVel * myEtaVolVel);
    //  New volume thermostat potential energy
    Real VEtaV = myEtaV * kbT;
    //  New box volume kinetic energy.
    Real pVol = (W * 0.5) * (myEpsilonVel * myEpsilonVel);
    //  New box volume potential energy.
    Real VVol = (myVolume * myTargetPres) / Constant::PRESSUREFACTOR;
    
    //  Add the energy from the extended system barostat.
    (app->energies)[ScalarStructure::INTEGRATOR] += pEtaV + VEtaV + pVol + VVol;
  }

  //  -----------------------------------------------------------------------  //
  //  doHalfKick().
  void NPTVerletIntegrator::doHalfKick() {

    //  Timestep.  Units: (fs)
    const Real halfDeltaT = 0.5 * getTimestep();

    // calculate the COM and momentum of each molecule
    buildMolecularCenterOfMass(&app->positions,app->topology);
    buildMolecularMomentum(&app->velocities,app->topology);

    // ---------------------------------------------------------------------
    //  Do the first update of the atom velocities
    //  Loop over all molecules
    for (unsigned int i=0; i < app->topology->molecules.size(); i++) {    
      //  Temporary storage element for the updated molecular momentum
      Vector3D Momentum(0.0,0.0,0.0);

      //  Loop over the atoms on this molecule
      for (unsigned int a=0; a < app->topology->molecules[i].size(); a++) {

        //  Current atom # and mass
        int atom = app->topology->molecules[i][a];
        Real mass = app->topology->atoms[atom].scaledMass; 

        //  Advance the velocities due to the thermostat and barostat forces.
        (app->velocities)[atom] *= exp(-(mass / app->topology->molecules[i].mass * (1.0 + 3.0/app->topology->degreesOfFreedom)
				       * myEpsilonVel + myEtaVel ) * halfDeltaT);

        //  Add to the new momentum of this molecule
        Momentum += (app->velocities)[atom] * mass;

      } //  end loop over atoms
      
      //  Store the updated molecular momentum
      app->topology->molecules[i].momentum = Momentum;
    } //  end loop over molecules

    // -------------------------------------------------------------------------
    //  Do the second and third updates of the atom velocities
    //  Loop over all molecules
    for (unsigned int i=0; i < app->topology->molecules.size(); i++) {
       
      //  Loop over the atoms on this molecule
      for (unsigned int a=0; a < app->topology->molecules[i].size(); a++) {
  
        //  Current atom # and mass
        int atom = app->topology->molecules[i][a];
        Real mass =  app->topology->atoms[atom].scaledMass;

        //  Advance the velocities due to barostat force.
        (app->velocities)[atom] -= (app->topology->molecules[i].momentum - (app->velocities)[atom] * mass)
	  * (1.0 + 3.0/app->topology->degreesOfFreedom) * (1.0 / app->topology->molecules[i].mass)
	  * myEpsilonVel * halfDeltaT;

        //  Advance the velocities due to atomic forces.
        (app->velocities)[atom] += (*myForces)[atom] * halfDeltaT * Constant::INV_TIMEFACTOR / mass;
      }  // end loop over atoms
    } //  end loop over molecules
  }  //  End doHalfKick().

  //  -----------------------------------------------------------------------  //
  //  doDrift().
  void NPTVerletIntegrator::doDrift() {

    //  Timestep.  Units: (fs)
    const Real deltaT = getTimestep();

    // -------------------------------------------------------------------------
    //  Do the first update of the atom positions
    //  Loop over all molecules
    for (unsigned int i=0; i < app->topology->molecules.size(); i++) {
       
      //  Loop over the atoms on this molecule
      for (unsigned int a=0; a < app->topology->molecules[i].size(); a++) {
  
        //  Current atom # and mass
        int atom = app->topology->molecules[i][a];
        Real mass = app->topology->atoms[atom].scaledMass;

        //  Advance the positions due to change in box length.
        (app->positions)[atom] *= exp((mass / app->topology->molecules[i].mass)
				    * myEpsilonVel * deltaT);
      }  // end loop over atoms
    } //  end loop over molecules

    // -------------------------------------------------------------------------
    //  Do the second and final update of the positions
    //  Loop over all molecules
    for (unsigned int i=0; i < app->topology->molecules.size(); i++) {

      //  Temporary storage element for the updated molecular COM
      Vector3D COM(0.0,0.0,0.0);

      //  Loop over the atoms on this molecule
      for (unsigned int a=0; a < app->topology->molecules[i].size(); a++) {
  
        //  Current atom # and mass
        int atom = app->topology->molecules[i][a];
        Real mass =  app->topology->atoms[atom].scaledMass;

        //  Advance the positions due to box volume or barostat (similar to COM shifting).
        (app->positions)[atom] += (app->topology->molecules[i].position - (app->positions)[atom] 
				 * (mass / app->topology->molecules[i].mass)) * deltaT * myEpsilonVel;

        //  Advance the positions due to velocity.
        (app->positions)[atom] += (app->velocities)[atom] * deltaT * Constant::INV_TIMEFACTOR;

        //  Add to the new COM of this molecule
        COM += (app->positions)[atom] * mass;
      }  //  end loop over atoms
      
      //  Store the updated molecular COM
      app->topology->molecules[i].position = COM / (app->topology->molecules[i].mass);
    } //  end loop over molecules

    //  update the COM of each molecule
    buildMolecularCenterOfMass(&app->positions,app->topology);
  }  // end doDrift()
  
  //  -------------------------------------------------------------------  //
  //  do2ndHalfkick()
  void NPTVerletIntegrator::do2ndHalfKick() {

    //  Timestep.  Units: (fs)
    const Real halfDeltaT = 0.5 * getTimestep();

    // -----------------------------------------------------------------------------
    //  Do the first update of the atomic velocities
    //  Loop over all molecules
    for (unsigned int i=0; i < app->topology->molecules.size(); i++) {

      //  Temporary storage element for the updated molecular momentum
      Vector3D Momentum(0.0,0.0,0.0);

      //  Loop over the atoms on this molecule
      for (unsigned int a=0; a < app->topology->molecules[i].size(); a++) {
  
        //  Current atom # and mass
        int atom = app->topology->molecules[i][a];
        Real mass =  app->topology->atoms[atom].scaledMass;

        //  Advance the velocities due to atomic forces.
        (app->velocities)[atom] += (*myForces)[atom] * halfDeltaT * Constant::INV_TIMEFACTOR / mass;

        //  Add to the new momentum of this molecule
        Momentum += (app->velocities)[atom] * mass;
      }  //  end loop over atoms

      //  Store the updated molecular momentum
      app->topology->molecules[i].momentum = Momentum;
    } //  end loop over molecules

    // ------------------------------------------------------------------------------
    //  Do the second and third updates of the atomic velocities
    //  Loop over all molecules
    for (unsigned int i=0; i < app->topology->molecules.size(); i++) {
       
      //  Temporary storage element for the updated molecular momentum
      Vector3D Momentum(0.0,0.0,0.0);

      //  Loop over the atoms on this molecule
      for (unsigned int a=0; a < app->topology->molecules[i].size(); a++) {
  
        //  Current atom # and mass
        int atom = app->topology->molecules[i][a];
        Real mass = app->topology->atoms[atom].scaledMass;

        //  Advance the velocities due to barostat force.
        (app->velocities)[atom] -= (app->topology->molecules[i].momentum - (app->velocities)[atom] * mass)
	  * (1.0 + 3.0/app->topology->degreesOfFreedom) * (1.0 / app->topology->molecules[i].mass)
	  * myEpsilonVel * halfDeltaT;

        //  Advance the velocities due to the thermostat and barostat forces.
        (app->velocities)[atom] *= exp(-(mass / app->topology->molecules[i].mass * (1.0 + 3.0/app->topology->degreesOfFreedom)
				       * myEpsilonVel + myEtaVel ) * halfDeltaT);

        //  Add to the new momentum of this molecule
        Momentum += (app->velocities)[atom] * mass;                                 
      }  // end loop over atoms

      //  Store the updated molecular momentum
      app->topology->molecules[i].momentum = Momentum;
    } //  end loop over molecules

    //  Add the new kinetic and potential energy of the system.
    Real PE = app->energies.potentialEnergy();
    Real KE = kineticEnergy(app->topology, &app->velocities);
    (app->energies)[ScalarStructure::INTEGRATOR] += KE + PE;
    
  }  // end do2ndHalfKick

  //  -----------------------------------------------------------------------  //
  //  run().
  void NPTVerletIntegrator::run(int numTimesteps) {
    for (int i = 0; i < numTimesteps; i++) {
      preStepModify();
      doHalfKick();
      doDriftOrNextIntegrator();
      calculateForces();
      do2ndHalfKick();
      postStepModify();
    }
  }
          
  //  -----------------------------------------------------------------------  //
  //  initialize().
  void NPTVerletIntegrator::initialize(ProtoMolApp* app) {

    STSIntegrator::initialize(app);

    buildMolecularCenterOfMass(&app->positions,app->topology);

    // initialize all forces and modifiers
    myForces->zero(app->positions.size());
    initializeForces();

    //  Initialize variables to something sane.
    myEpsilonVel = 0.;
    myEta        = 0.;
    myEtaV       = 0.;
    myEtaVel     = 0.;
    myEtaVolVel  = 0.;

    // get the initial box volume and # of atoms
    myVolume = app->topology->getVolume(app->positions);

    //  Compute the fixed barostat mass and thermostat masses.
    Qo = app->topology->degreesOfFreedom * kbT * (myTauT * myTauT);
    Qv = kbT * (myTauV * myTauV);
    W = (app->topology->degreesOfFreedom + 3) * kbT * (myTauP * myTauP);

  }

  //  -----------------------------------------------------------------------  //
  //  createRattleModifier().
  Modifier* NPTVerletIntegrator::createRattleModifier(Real eps, int maxIter){
    return (new ModifierNPTRattle<NPTVerletIntegrator>(eps, maxIter, this));
  }
  
  //  -----------------------------------------------------------------------  //
  //  createShakeModifier().
  Modifier* NPTVerletIntegrator::createShakeModifier(Real eps, int maxIter){
    return (new ModifierNPTShake<NPTVerletIntegrator>(eps, maxIter, this));
  }

  //  -----------------------------------------------------------------------  //
  //  Add the modifiers
  void NPTVerletIntegrator::addModifierAfterInitialize() {
    adoptPreStepModifier(new ModifierPreForceThermostat<NPTVerletIntegrator>(this,1)); 
    adoptPreStepModifier(new ModifierPreForceBarostat<NPTVerletIntegrator>(this,2));
    adoptPostStepModifier(new ModifierPostForceBarostat<NPTVerletIntegrator>(this,1));
    adoptPostStepModifier(new ModifierPostForceThermostat<NPTVerletIntegrator>(this,2));
    STSIntegrator::addModifierAfterInitialize();
  }

  //  -----------------------------------------------------------------------  //
  //  getParameters().
  void NPTVerletIntegrator::getParameters(vector<Parameter>& parameters) const {
    STSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("temperature", Value(myTargetTemp)));
    parameters.push_back(Parameter("pressure", Value(myTargetPres)));
    parameters.push_back(Parameter("tauT", Value(myTauP)));
    parameters.push_back(Parameter("tauV", Value(myTauV)));
    parameters.push_back(Parameter("tauP", Value(myTauP)));
  }

  //  -----------------------------------------------------------------------  //
  //  doMake().
  STSIntegrator* NPTVerletIntegrator::doMake(const vector<Value>& values,ForceGroup* fg)const{
    return new NPTVerletIntegrator(values[0],values[1],values[2],values[3],values[4],values[5],fg);
  }

}
