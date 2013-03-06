/* -*- c++ -*- */
#ifndef FORCE_FACTORY_H
#define FORCE_FACTORY_H

#include <protomol/base/Factory.h>
#include <protomol/force/Force.h>
#include <protomol/config/Value.h>

namespace ProtoMol {
  //________________________________________ ForceFactory
  class ForceFactory : public Factory<Force> {
    typedef std::set<std::string, ltstrNocase> policy_t;
    typedef std::map<std::string, std::string, ltstrNocase> policiesSorted_t;

    struct ForceType {
      policy_t policy;
      policy_t policies;
      policiesSorted_t policiesSorted;
    };

    typedef std::map<std::string, ForceType, ltstrNocase> forceTypes_t;
    mutable forceTypes_t forceTypes;
    typedef std::map<std::string, std::string, ltstrNocase>
    forceTypesSorted_t;
    mutable forceTypesSorted_t forceTypesSorted;
    mutable CompareForce *lastCompareForce;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ForceFactory();
    virtual ~ForceFactory();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From FactoryBase<Force>
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    virtual std::ostream &print(std::ostream &stream) const;
    virtual void registerHelpText() const {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class DetailsForceFactory
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Force *make(const std::string &id, std::vector<Value> values =
                std::vector<Value>()) const;

  private:
    void updateCache() const;
    std::vector<std::string> splitForceString(const std::string &id) const;
    std::vector<std::string>
    splitForceStringSorted(const std::string &id) const;
    std::string sortForceString(const std::string &id) const;
    std::string uniqueForceString(const std::string &id) const;
  };
}
#endif /* FORCE_FACTORY_H */
