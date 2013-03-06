/*  -*- c++ -*-  */
#ifndef MODIFIER_H
#define MODIFIER_H
#include <vector>

#include <protomol/base/Makeable.h>

#include <ostream>

namespace ProtoMol {
  class ProtoMolApp;
  class Vector3DBlock;
  class Integrator;

  //________________________________________ Modifier
  /**
     Base class for all kind of modifier implementation, based on the Observer
     Pattern. A modifier object can be added at five different stages of a
     single time step:@n
     - before doDriftOrNextIntegrator (pre) @n
     - after doDriftOrNextIntegrator (post) @n
     - before the force calculation (pre) @n
     - between system and extended force calculation (medi) @n
     - after the force calculation (post) @n @n

     Each modifier object performs it changes to the system via the method 
     execute(...) with the implementation details in doExecute(). The execution
     is defined 1. by an order number and 2. by the pointer of the object.
     An implementation can be either internal -- added by an integrator -- or
     external -- added by the user or at application level. Integrator provides
     all methods add, remove and delete modifiers. Integrators provides also
     interface to add internal modifiers before/after the initialization of
     forces.

     Internal modifiers (isInternal() is true) are removed, if any, by
     Intergrator under initialize and each integrator has to add its modifiers,
     if any. To add internal modifier you should override
     addModifierBeforeInitialize() and/or addModifierAfterInitialize() in order
     the specify if the modifications should be considered during the (force)
     initialization or not. Furthermore, it is possible to disable and enable
     a modifier.
   */
  class Modifier : public Makeable<Modifier> {
  public:
    static std::string scope;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Modifier(int order = 0) : 
      myOrder(order), myEnable(true), app(0), myForces(0) {}
    virtual ~Modifier() {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getScope() const {return scope;}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class Modifier
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    /// The method, which calls the implemenation
    void execute(Integrator *i);

    /// If the modifier is internal (added by an integrator) or
    /// external (added by the user)
    virtual bool isInternal() const = 0;

    /// Returns order of execution
    int order() const {return myOrder;}

    /// Activate modifier
    void enable() const {myEnable = true;}
    /// Deactivate modifier
    void disable() const {myEnable = false;}
    /// If the modifier is active (doExecute() is called)
    bool isEnabled() const {return myEnable;}

    /// Strict weak order using first order and than pointer to use set<>
    bool operator<(const Modifier &m) const {
      if (myOrder < m.myOrder) return true;
      else if (myOrder > m.myOrder) return false;
      return this < &m;
    }

    /// Initialize
    void initialize(ProtoMolApp *app, Vector3DBlock *forces);

    virtual Modifier *doMake(const std::vector<Value> &values) const {return 0;}
    virtual void getParameters(std::vector<Parameter> &parameters) const {}

    /// print/debug
    virtual std::ostream &print(std::ostream &stream) const {
      return stream << getId() << endl;
    }

    friend std::ostream &operator<<(std::ostream &, const Modifier &);

  private:
    /// The method, which does the actual modification
    virtual void doExecute(Integrator *i) = 0;
    /// Implemenation of initialize
    virtual void doInitialize() {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    int myOrder;
    mutable bool myEnable;

  protected:
    ProtoMolApp *app;
    Vector3DBlock *myForces;
  };

  inline std::ostream &operator<<(std::ostream &stream, const Modifier &m) {
    return m.print(stream);
  }
}
#endif /* MODIFIER_H */
