/*  -*- c++ -*-  */
#ifndef MAKABLEDEFINITION_H
#define MAKABLEDEFINITION_H

#include <vector>
#include <protomol/config/Parameter.h>

namespace ProtoMol {
  //________________________________________ MakeableDefinition
  /**
   * Vector container struct for object definitions.
   */
  struct MakeableDefinition {
    // Container struct for makable definitions

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    MakeableDefinition() {}
    MakeableDefinition(const std::string &i, const std::vector<Parameter> &p) :
      id(i), parameters(p) {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    std::string id;         ///< keyword of the object
    std::vector<Parameter> parameters; ///< parameters of the obejct
  };
}
#endif /* MAKABLEDEFINITION_H */
