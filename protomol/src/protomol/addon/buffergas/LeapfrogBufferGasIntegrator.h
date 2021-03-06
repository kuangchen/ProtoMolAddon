/*  -*- c++ -*-  */
#ifndef __LEAPFROG_BUFFERGAS_INTEGRATOR_H_
#define __LEAPFROG_BUFFERGAS_INTEGRATOR_H_

#include <protomol/integrator/STSIntegrator.h>
#include <string>
#include <vector>
                
namespace ProtoMolAddon {
  namespace BufferGas {

    using namespace ProtoMol;

    template <class Collision>
    class LeapfrogBufferGasIntegrator : public STSIntegrator {

    private:
      std::string fname;
      Collision collision;

    public:
      LeapfrogBufferGasIntegrator(): 
	STSIntegrator(), 
	collision()
	{}
      
      LeapfrogBufferGasIntegrator(Real timestep, const std::string &fname, 
				  ForceGroup *overloadedForces) :
	STSIntegrator(timestep, overloadedForces),
	collision(fname)
	{}

      ~LeapfrogBufferGasIntegrator() {}
      
    public:
      virtual void initialize(ProtoMolApp *app) {
	STSIntegrator::initialize(app);
	collision.Initialize(app);
	initializeForces();
      }

      void doHalfKickdoDrift() {
	if (anyPreDriftOrNextModify()) {
	  doHalfKick();
	  doDriftOrNextIntegrator();
	} else {
	  Real h = getTimestep() * INV_TIMEFACTOR;
	  const unsigned int count = app->positions.size();

	  //  Do a half kick on beta.
	  updateBeta(0.5 * h);

	  for (unsigned int i = 0; i < count; ++i) {
	    app->velocities[i] += (*myForces)[i] * h * 0.5 /
	      app->topology->atoms[i].scaledMass;
	    // app->positions[i] += app->velocities[i] * h;
	  }
    
	  app->positions += app->velocities*h;

	  buildMolecularCenterOfMass(&app->positions, app->topology);
	  buildMolecularMomentum(&app->velocities, app->topology);
	  postDriftOrNextModify();
	}
      }

      void doKickdoDrift() {
	if (anyPreDriftOrNextModify() || anyPreStepModify() ||
	    anyPostStepModify()) {
	  if (anyPreStepModify() || anyPostStepModify()) {
	    doHalfKick();
	    postStepModify();
	    preStepModify();
	    doHalfKick();
	  } else
	    doKick();
	  doDriftOrNextIntegrator();
	} else {
	  Real h = getTimestep() * INV_TIMEFACTOR;
	  const unsigned int count = app->positions.size();

	  updateBeta(h);

	  for (unsigned int i = 0; i < count; ++i) {
	    app->velocities[i] +=
	      (*myForces)[i] * h / app->topology->atoms[i].scaledMass;
	    //  app->positions[i] += app->velocities[i] * h;
	  }

	  app->positions += app->velocities*h;

	  buildMolecularCenterOfMass(&app->positions, app->topology);
	  buildMolecularMomentum(&app->velocities, app->topology);
	  postDriftOrNextModify();
	}
      }

      void run(int numTimesteps) {
	double dt = getTimestep() / SI::TIME_FS;
	if (numTimesteps < 1)
	  return;

	preStepModify();

	collision.Collide(dt);

	doHalfKickdoDrift();
	calculateForces();

	for (int i = 1; i < numTimesteps; i++) {
	  collision.Collide(dt);
	  doKickdoDrift();
	  calculateForces();
	}

	doHalfKick();
	postStepModify();
      }

      STSIntegrator* doMake(const std::vector<Value> &values, ForceGroup *fg) const {
	return new LeapfrogBufferGasIntegrator<Collision>(values[0], values[1], fg);
      }

//  --------------------------------------------------------------------  //
//  This function is necessary to compute the shadow Hamiltonian and it   //
//  is integrator specific.  This version is written to work with LF.     //
//  Update beta: beta -= dt * ( q * F + 2 U )                             //
//  --------------------------------------------------------------------  //

      void updateBeta(Real dt) {
	//  ----------------------------------------------------------------  //
	//  The shadow calculation is done in a postStep modifier.  If there  //
	//  aren't any, then obviously we don't need to do this calculation.  //
	//  It's possible that a different poststep modifier could make this  //
	//  execute, but no harm would be done ... only some extra cycles.    //
	//  ----------------------------------------------------------------  //

	if (!(anyPostStepModify() || top()->anyPostStepModify()))
	  return;

	Real posDotF = 0.;

	for (unsigned int i = 0; i < app->positions.size(); i++)
	  posDotF += app->positions[i].dot((*myForces)[i]);

	myBeta -= dt * (posDotF + 2. * myPotEnergy);
      }

      void getParameters(std::vector<Parameter> &parameters) const {
	STSIntegrator::getParameters(parameters);
	parameters.push_back(Parameter("filename", Value(fname, ConstraintValueType::NotEmpty())));
      }

      virtual std::string getIdNoAlias() const 
	{return std::string("LeapfrogBufferGas") + Collision::GetKeyword();}
    };

  }
}

#endif

