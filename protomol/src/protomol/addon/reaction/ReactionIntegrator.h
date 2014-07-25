#ifndef __REACTION_INTEGRATOR_H
#define __REACTION_INTEGRATOR_H

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

#include <protomol/addon/reaction/Reaction.h>
#include <protomol/integrator/STSIntegrator.h>


namespace ProtoMolAddon {
  namespace Reaction {

    using namespace ProtoMol;
    using namespace ProtoMol::Report;
    
    template <class Mechanism>
    class ReactionIntegrator : public STSIntegrator {

    public:
      ReactionIntegrator() : STSIntegrator(), r() {}
	
      ReactionIntegrator(Real timestep, const std::string &fname, ForceGroup *overloadedForces) : STSIntegrator(timestep, overloadedForces), r(fname) {}
   

    protected:
      void doKickdoDrift();
      void doHalfKickdoDrift();

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // From class Makeable
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
      virtual std::string getIdNoAlias() const {return Mechanism::keyword+"Integrator";}
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
      Reaction<Mechanism> r;

    public:
      static const std::string keyword;

    };

    template<class Mechanism>
    const std::string ReactionIntegrator<Mechanism>::keyword("Integrator");

    template<class Mechanism>
    void ReactionIntegrator<Mechanism>::initialize(ProtoMolApp *app) {
      STSIntegrator::initialize(app);
      initializeForces();
      r.init(app);
    }

    template<class Mechanism>
    void ReactionIntegrator<Mechanism>::doHalfKickdoDrift() {
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

	r.react(getTimestep() * Constant::ToSI::time * 0.5);
	buildMolecularCenterOfMass(&app->positions, app->topology);
	buildMolecularMomentum(&app->velocities, app->topology);
	postDriftOrNextModify();
      }
    }

    template<class Mechanism>
    void ReactionIntegrator<Mechanism>::doKickdoDrift() {
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

	r.react(getTimestep() * Constant::ToSI::time);
	buildMolecularCenterOfMass(&app->positions, app->topology);
	buildMolecularMomentum(&app->velocities, app->topology);
	postDriftOrNextModify();
      }
    }

    template<class Mechanism>
    void ReactionIntegrator<Mechanism>::run(int numTimesteps) {
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

    template<class Mechanism>
    STSIntegrator *ReactionIntegrator<Mechanism>::doMake(const vector<Value> &values, ForceGroup *fg) const {
      return new ReactionIntegrator(values[0], values[1], fg);
    }

//  --------------------------------------------------------------------  //
//  This function is necessary to compute the shadow Hamiltonian and it   //
//  is integrator specific.  This version is written to work with LF.     //
//  Update beta: beta -= dt * ( q * F + 2 U )                             //
//  --------------------------------------------------------------------  //
    template<class Mechanism>
    void ReactionIntegrator<Mechanism>::updateBeta(Real dt) {
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

    template<class Mechanism>
    void ReactionIntegrator<Mechanism>::getParameters(vector<Parameter> &parameters) const {
      STSIntegrator::getParameters(parameters);
      parameters.push_back(Parameter("filename", Value(filename, ConstraintValueType::NotEmpty())));
    }


  }
}

#endif
