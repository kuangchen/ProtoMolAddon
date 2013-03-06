/*  -*- c++ -*-  */
#ifndef MODIFIERPREFORCEBAROSTAT_H
#define MODIFIERPREFORCEBAROSTAT_H

#include <protomol/modifier/Modifier.h>

namespace ProtoMol {

  //_________________________________________________________________ ModifierPreForceBarostat
  template<class TIntegrator>
  class ModifierPreForceBarostat : public Modifier {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ModifierPreForceBarostat(TIntegrator* i):Modifier(Constant::MAX_INT-400),myTheIntegrator(i){}
    ModifierPreForceBarostat(TIntegrator* i, int order):Modifier(order),myTheIntegrator(i){}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Modifier
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool isInternal() const {return true;}
  private:
    virtual void doExecute(Integrator* myTheIntegrator){
      ((TIntegrator*)myTheIntegrator)->PreForceBarostat();
    }
    virtual std::string getIdNoAlias()const{return std::string("PreForceBarostat");};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    TIntegrator* myTheIntegrator;
  };

}
#endif /* MODIFIERPREFORCEBAROSTAT_H */
