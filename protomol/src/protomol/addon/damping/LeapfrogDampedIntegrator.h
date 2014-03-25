/*  -*- c++ -*-  */
#ifndef LEAPFROGDAMPEDINTEGRATOR_H
#define LEAPFROGDAMPEDINTEGRATOR_H

#include <string>
#include <protomol/integrator/STSIntegrator.h>
#include <protomol/addon/damping/DampingSpec.h>

#include <set>
#include <vector>

using ProtoMol::Integrator;
using ProtoMolAddon::Damping::DampingSpec;
using std::vector;
using std::set;
using std::string;
using namespace ProtoMol;

namespace ProtoMolAddon {
  namespace Damping {
  
    //____ LeapfrogDampedIntegrator
    class LeapfrogDampedIntegrator : public STSIntegrator {
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Constructors, destructors, assignment
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
      LeapfrogDampedIntegrator();
      LeapfrogDampedIntegrator(Real timestep, const string &spec_fname,
			       ForceGroup *overloadedForces);
      ~LeapfrogDampedIntegrator();

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // New methods of class LeapfrogDampedIntegrator
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    protected:
      void doKickdoDrift();
      void doHalfKickdoDrift();

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // From class Makeable
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
      virtual std::string getIdNoAlias() const {return keyword;}
      virtual void getParameters(std::vector<Parameter> &parameters) const ;
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
    public:
      static const std::string keyword;

    private:
      string spec_fname;
      DampingSpec spec;
    };

  }
}

#endif

