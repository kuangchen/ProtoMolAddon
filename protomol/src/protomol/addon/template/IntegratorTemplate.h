#ifndef _INTEGRATOR_TEMPLATE_H
#define _INTEGRATOR_TEMPLATE_H

#include <vector>
#include <string>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/module/MainModule.h>
#include <protomol/integrator/STSIntegrator.h>

namespace ProtoMolAddon {
  namespace Template {

    using namespace ProtoMol;
    using namespace ProtoMol::Report;
    
    template <class Reaction>
    class GenericIntegrator : public STSIntegrator {

    public:
      GenericIntegrator() : STSIntegrator(), r() {}
      GenericIntegrator(Real timestep, const std::string &fname, ForceGroup *overloadedForces) : STSIntegrator(timestep, overloadedForces), r(fname) {}

    protected:
      void doKickdoDrift();
      void doHalfKickdoDrift();

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // From class Makeable
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
      virtual std::string getIdNoAlias() const {return Reaction::GetName(); }
      virtual void getParameters(std::vector<Parameter> &parameters) const;
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // From class Integrator
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
      virtual void initialize(ProtoMolApp *app);
      virtual void run(int numTimesteps);
      virtual void updateBeta(Real dt);
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // From class STSIntegrator
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    private:
      virtual STSIntegrator *doMake(const std::vector<Value> &values,
				    ForceGroup *fg) const;
      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // My data members
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    private:
      std::string filename; 
      Reaction r;

    public:
      static const std::string keyword;

    };

    template<class Reaction>
    const std::string GenericIntegrator<Reaction>::keyword("Integrator");

    template<class Reaction>
    void GenericIntegrator<Reaction>::initialize(ProtoMolApp *app) {
      STSIntegrator::initialize(app);
      initializeForces();
      r.Initialize(app);
    }

    template<class Reaction>
    void GenericIntegrator<Reaction>::doHalfKickdoDrift() {
      if (anyPreDriftOrNextModify()) {
	doHalfKick();
	doDriftOrNextIntegrator();
      } else {
	Real h = getTimestep() * ProtoMol::Constant::INV_TIMEFACTOR;
	const unsigned int count = app->positions.size();

	//  Do a half kick on beta.
	updateBeta(0.5 * h);

	for (unsigned int i = 0; i < count; ++i) {
	  app->velocities[i] += (*myForces)[i] * h * 0.5 /
	    app->topology->atoms[i].scaledMass;
	  // app->positions[i] += app->velocities[i] * h;
	}
    
	app->positions += app->velocities*h;

	r.Update(app->topology->time * Constant::ToSI::time,
		 getTimestep() * Constant::ToSI::time * 0.5);

	buildMolecularCenterOfMass(&app->positions, app->topology);
	buildMolecularMomentum(&app->velocities, app->topology);
	postDriftOrNextModify();
      }
    }

    template<class Reaction>
    void GenericIntegrator<Reaction>::doKickdoDrift() {
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
	Real h = getTimestep() * ProtoMol::Constant::INV_TIMEFACTOR;

	const unsigned int count = app->positions.size();

	updateBeta(h);

	for (unsigned int i = 0; i < count; ++i) {
	  app->velocities[i] +=
	    (*myForces)[i] * h / app->topology->atoms[i].scaledMass;
	  //  app->positions[i] += app->velocities[i] * h;
	}

	app->positions += app->velocities*h;

	r.Update(app->topology->time * Constant::ToSI::time,
		 getTimestep() * Constant::ToSI::time);
	
	buildMolecularCenterOfMass(&app->positions, app->topology);
	buildMolecularMomentum(&app->velocities, app->topology);
	postDriftOrNextModify();
      }
    }

    template<class Reaction>
    void GenericIntegrator<Reaction>::run(int numTimesteps) {
      if (numTimesteps < 1)
	return;

      preStepModify();
      doHalfKickdoDrift();
      calculateForces();

      for (int i = 1; i < numTimesteps; i++) {
	doKickdoDrift();
	calculateForces();

      }
      doHalfKick();
      postStepModify();
    }

    template<class Reaction>
    STSIntegrator *GenericIntegrator<Reaction>::doMake(const vector<Value> &values, ForceGroup *fg) const {
      return new GenericIntegrator(values[0], values[1], fg);
    }

//  --------------------------------------------------------------------  //
//  This function is necessary to compute the shadow Hamiltonian and it   //
//  is integrator specific.  This version is written to work with LF.     //
//  Update beta: beta -= dt * ( q * F + 2 U )                             //
//  --------------------------------------------------------------------  //
    template<class Reaction>
    void GenericIntegrator<Reaction>::updateBeta(Real dt) {
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

    template<class Reaction>
    void GenericIntegrator<Reaction>::getParameters(vector<Parameter> &parameters) const {
      STSIntegrator::getParameters(parameters);
      parameters.push_back(Parameter(r.GetParameterName(), Value(filename, ConstraintValueType::NotEmpty())));
    }


  }
}

#endif
