#ifndef __IONIZATION_INTEGRATOR_H_
#define __IONIZATION_INTEGRATOR_H_

#include <vector>
#include <string>
#include <protomol/integrator/STSIntegrator.h>
#include <protomol/addon/ionization/IonizationManager.h>

using namespace ProtoMol;

namespace ProtoMolAddon {
  namespace Ionization {
    
    
    class IonizationIntegrator : public STSIntegrator {
      
    public:
      IonizationIntegrator();
      IonizationIntegrator(Real timestep, const std::string &fname, ForceGroup *overloadedForces);
      

    protected:
      void doKickdoDrift();
      void doHalfKickdoDrift();

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // From class Makeable
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
      virtual std::string getIdNoAlias() const {return keyword;}
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
      IonizationManager im;

    public:
      static const std::string keyword;
      
    };
  }
}

#endif
