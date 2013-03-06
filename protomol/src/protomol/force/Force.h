/* -*- c++ -*- */
#ifndef FORCE_H
#define FORCE_H

#include <protomol/base/Makeable.h>
#include <protomol/base/MakeableDefinition.h>
#include <protomol/base/Exception.h>

#include <vector>

namespace ProtoMol {
  class ForceGroup;
  class ScalarStructure;
  class GenericTopology;
  class Vector3DBlock;
  class CompareForce;
  class TimeForce;
  class ForceGroup;
  //________________________________________ Force

  class Force : public Makeable<Force> {
    // This class contains the definition of one force

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Force() {}
    virtual ~Force() {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class Force
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    virtual unsigned int numberOfBlocks(const GenericTopology *,
                                        const Vector3DBlock *) {return 1;}
    //< Returns the number of blocks the actual force can be split up. To
    //< run it sequential return 1. The number of blocks must match with the
    //< number of next() calls in parallelEvaluate(){...} running with > 1 node.

    virtual std::string getKeyword() const = 0;
    //< Returns the actual keyword of the concrete force:
    //<
    //< Your force *Force has private inheritance from *ForceBase
    //< and needs following virtual method:
    //< virtual std::string getKeyword() const{return keyword;}
    //<
    //< Define all static and basics, which do not depend on templates in
    //< *ForceBase.

    Force *make(const std::vector<Value> &values) const;
    //< Factory method

    virtual void uncache() {};
    //< Marks the cache as out of date, forces to call initialize

    virtual void addToForceGroup(ForceGroup *forceGroup) = 0;
    //< Adds the force to corresponding force subgroups System or Extended

    virtual CompareForce *makeCompareForce(Force *actualForce,
                                           CompareForce *compareForce) const;
    //< Creating the right instance of a ComparForce object.

    virtual TimeForce *makeTimeForce(Force *actualForce) const;
    //< Creating the right instance of a TimeForce object.

    void setParameters(std::vector<Value> value);
    //< update of parameters

    template<class T>
    void setParameter(const std::string &key, T val) {
      std::vector<Parameter> parameters;
      getParameters(parameters);
      std::vector<Value> values(parameters.size());
      bool update = false;
      for (unsigned int i = 0; i < parameters.size(); i++) {
        values[i] = parameters[i].value;
        if (equalNocase(key,
              parameters[i].keyword) && parameters[i].value != val) {
          values[i] = val;
          update = true;
        }
      }

      if (update) setParameters(values);
    }
    //< update of parameters with given value which match the keyword

    template<class T>
    void setParameter(unsigned int index, T val) {
      std::vector<Parameter> parameters;
      getParameters(parameters);
      std::vector<Value> values(parameters.size());
      bool update = false;
      for (unsigned int i = 0; i < parameters.size(); i++) {
        values[i] = parameters[i].value;
        if (i == index && parameters[i].value != val) {
          values[i] = val;
          update = true;
        }
      }

      if (update) setParameters(values);
    }
    //< update of parameter index with given value

  private:
    virtual void doSetParameters(std::vector<Value>) {}

  protected:
    template<class T>
    void setParametersBySwapping(T *obj, std::vector<Value> values) {
      T *tmp = dynamic_cast<T *>(obj->make(values));
      if (tmp != NULL) {
        std::swap(*obj, *tmp);
        delete tmp;
      }
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getScope() const {return scope;}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string scope;
  };

  //________________________________________ INLINES
}
#endif /* FORCE_H */
