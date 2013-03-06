#ifndef PYMTSINTEGRATOR_H
#define PYMTSINTEGRATOR_H


#include <protomol/integrator/MTSIntegrator.h>

namespace ProtoMol {

class PyMTSIntegrator : public MTSIntegrator {
 public:
  PyMTSIntegrator() : MTSIntegrator() {}
  PyMTSIntegrator(int cycles, ForceGroup *overloadedForces,
                  StandardIntegrator *nextIntegrator) : 
    MTSIntegrator(cycles, overloadedForces, nextIntegrator) {}
 private:
    virtual MTSIntegrator *doMake(const std::vector<Value> &values,
				    ForceGroup *fg, StandardIntegrator *nextIntegrator) const {return NULL;}


};

}




#endif
