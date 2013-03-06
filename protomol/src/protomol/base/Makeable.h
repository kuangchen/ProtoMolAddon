/*  -*- c++ -*-  */
#ifndef PROTOMOL_MAKEABLE_H
#define PROTOMOL_MAKEABLE_H

#include <protomol/base/MakeableDefinition.h>

#include <iostream>
using namespace std;

namespace ProtoMol {
  class Configuration;

  /***
      Base class of all object, which can be create dynamically based on a 
      prototype, normally used together with a Factory.
   */
  class MakeableBase {
    std::string alias;

  public:
    virtual ~MakeableBase() {};

    /// retrieve all parameters
    virtual void getParameters(std::vector<Parameter> &parameters) const = 0;
    virtual MakeableDefinition getDefinition() const;

    std::string getId() const;
    std::string getAlias() const;
    std::string setAlias(const std::string &id);
    virtual std::string getIdNoAlias() const = 0;
    virtual std::string getScope() const = 0;
    virtual std::string getText() const {return std::string();}

    void assertParameters(const std::vector<Value> &values) const;
    bool checkParameters(const std::vector<Value> &values) const;
    bool checkParameterTypes(const std::vector<Value> &values) const;
    virtual bool adjustWithDefaultParameters(std::vector<Value> &values,
                                             const Configuration *) const
    {return checkParameterTypes(values);}

    template<typename T>
    static T *copy(T *obj) {
      T *clone = NULL;
      if (obj != NULL) {
        std::string err;
        std::vector<Parameter> p;
        obj->getParameters(p);
        std::vector<Value> v(p.size());
        for (unsigned int i = 0; i < p.size(); i++)
          v[i] = p[i].value;

        clone = obj->make(err, v);
      }
      return clone;
    }

  protected:
    template<typename T>
    T *adjustAlias(T *obj) const {
      if (obj) obj->setAlias(getId());
      return obj;
    }
  };

  
  template <class T>
  class Makeable : public MakeableBase {
  public:
    virtual ~Makeable() {}

  protected:
    virtual T *doMake(const std::vector<Value> &values) const = 0;
  };
}
#endif // PROTOMOL_MAKEABLE_H
