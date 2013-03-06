/*  -*- c++ -*-  */
#ifndef MODIFIERNPTRATTLE_H
#define MODIFIERNPTRATTLE_H

#include <protomol/modifier/ModifierMetaRattle.h>
#include <protomol/type/Vector3DBlock.h>

namespace ProtoMol {

  //_________________________________________________________________ModifierNPTRattleDetails
  class ModifierNPTRattleDetails : public ModifierMetaRattle {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ModifierNPTRattleDetails(Real eps, int maxIter, int order);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class ModifierNPTRattleDetails
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual Real getEpsilonVel() const=0;
    virtual Real getEtaVel() const=0;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Modifier
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual void doExecute(Integrator*);
    virtual std::string getIdNoAlias() const {return std::string("NPTRattle");};
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
  };


  //_________________________________________________________________ModifierNPTRattle
  template<class TIntegrator>
  class ModifierNPTRattle : public ModifierNPTRattleDetails {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ModifierNPTRattle(Real eps, 
		      int maxIter,
		      const TIntegrator* i,
		      int order=Constant::MAX_INT-400):ModifierNPTRattleDetails(eps,maxIter,order),
						       myTheIntegrator(i){}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class ModifierMetaRattleShake
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    virtual Real getTimestep()const{return myTheIntegrator->getTimestep();}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class ModifierNPTRattleDetails
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual Real getEpsilonVel() const{return myTheIntegrator->getEpsilonVel();}
    virtual Real getEtaVel() const{return myTheIntegrator->getEtaVel();}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    const TIntegrator* myTheIntegrator;
  };

}
#endif /* MODIFIERNPTRATTLE_H */
