/*  -*- c++ -*-  */
#ifndef HERMITE_H
#define HERMITE_H

#include <protomol/type/Real.h>
#include <string>

namespace ProtoMol {
  //_________________________________________________________________ Hermite
  class Hermite {
    // Hermite interpolation.
    // theta are the lk, where dTheta are the derivatives of lk


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Hermite();
    Hermite(unsigned int order);
    Hermite(unsigned int order, Real w);
    ~Hermite();
    Hermite(const Hermite& Hermite);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class Hermite
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    unsigned int getOrder() const {return myInterOrder;}
    void setOrder(unsigned int order);
    void set(Real w);

    static bool isSigma(unsigned int order){return (Hermite(order,0.0).theta[order-1] == 0.0 && Hermite(order,0.0).dTheta[order-1] == 0.0);}
    // true iff one theta is 1 and all other 0 for w=0

    static const std::string& getKeyword() {return keyword;}
    // Returns the keyword/name of this interpolation

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    unsigned int myInterOrder;
  public:
    Real* theta;
    Real* dTheta;
  public:
    static const std::string keyword;
  };
  //______________________________________________________________________ INLINES

}
#endif /* HERMITE_H */
