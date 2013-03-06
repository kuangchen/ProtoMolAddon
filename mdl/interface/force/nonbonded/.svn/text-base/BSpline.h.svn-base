/*  -*- c++ -*-  */
#ifndef BSPLINE_H
#define BSPLINE_H

#include <protomol/type/Real.h>
#include <string>

namespace ProtoMol {
  //_________________________________________________________________ BSpline
  /**
   * BSpline interpolation.
   * theta are the Bk, where dTheta are the derivatives of Bk
   */
  class BSpline {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    BSpline();
    BSpline(unsigned int order);
    BSpline(unsigned int order, Real w);
    ~BSpline();
    BSpline(const BSpline& bspline);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class BSpline
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    unsigned int getOrder() const {return myInterOrder;}
    void setOrder(unsigned int order);
    void set(Real w);

    static bool isSigma(unsigned int order){return (BSpline(order,0.0).theta[order-1] == 0.0 && BSpline(order,0.0).dTheta[order-1] == 0.0);}
    ///< true iff one theta is 1 and all other 0 for w=0

    static const std::string& getKeyword() {return keyword;}
    ///< Returns the keyword/name of this interpolation

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    unsigned int myInterOrder;
  public:
    Real* theta; ///< Bk
    Real* dTheta;///< derivatives of Bk
  public:
    static const std::string keyword;
  };
  //______________________________________________________________________ INLINES
}
#endif /* BSPLINE_H */
