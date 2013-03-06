/* -*- c++ -*- */
#ifndef METAFORCE_H
#define METAFORCE_H

#include <protomol/force/Force.h>

namespace ProtoMol {
  class ForceGroup;
  //________________________________________ MetaForce

  class MetaForce : public Force {
    // This class contains the definition of one Meta force (a force that
    // works only on positions)

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    MetaForce() : Force() {}
    virtual ~MetaForce() {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class MetaForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void addToForceGroup(ForceGroup *forceGroup);
    virtual void getDeepForces(std::vector<Force *> &forces) const = 0;
  };

  //________________________________________ INLINES
}
#endif /* METAFORCE_H */



