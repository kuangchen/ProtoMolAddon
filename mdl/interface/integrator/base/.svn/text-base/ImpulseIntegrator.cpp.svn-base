#include "ImpulseIntegrator.h"
#include <protomol/type/Vector3DBlock.h>

using std::vector;
using std::string;

namespace ProtoMol {
  //_________________________________________________________________ ImpulseIntegrator
  const string ImpulseIntegrator::keyword("Impulse");

  ImpulseIntegrator::ImpulseIntegrator():MTSIntegrator(){}

  ImpulseIntegrator::ImpulseIntegrator(int cycles, 
				       ForceGroup *overloadedForces,
				       StandardIntegrator *nextIntegrator) :
    MTSIntegrator(cycles, overloadedForces, nextIntegrator) {

  }

  void ImpulseIntegrator::initialize(ProtoMolApp* app) {
    MTSIntegrator::initialize(app);
    initializeForces();
  }

  MTSIntegrator* ImpulseIntegrator::doMake(const vector<Value>& values,ForceGroup* fg, StandardIntegrator *nextIntegrator)const{
    return new ImpulseIntegrator(values[0],fg,nextIntegrator);
  }


    //  --------------------------------------------------------------------  //
    //  This function is necessary to compute the shadow Hamiltonian and it   //
    //  is integrator specific.  This version is written to work with imp.    //
    //  beta -= dt * ( q * (F_slow + F_fast) + 2 ( U_slow + U_fast ) )        //
    //  --------------------------------------------------------------------  //

    void ImpulseIntegrator::updateBeta( Real dt ) {

        //  ----------------------------------------------------------------  //
        //  The shadow calculation is done in a postStep modifier.  If there  //
        //  aren't any, then obviously we don't need to do this calculation.  //
        //  It's possible that a different poststep modifier could make this  //
        //  execute, but no harm would be done ... only some extra cycles.    //
        //  ----------------------------------------------------------------  //

        if( ! ( anyPostStepModify() || top()->anyPostStepModify() ) )
            return;

        Real posDotF = 0.;

        for( unsigned int i = 0; i < app->positions.size(); i++ )
            posDotF += (app->positions)[i].dot( (*myForces)[i] );

        myBeta -= dt * ( posDotF + 2. * myPotEnergy );

    }


}

