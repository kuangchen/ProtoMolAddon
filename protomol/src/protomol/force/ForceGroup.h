/* -*- c++ -*- */
#ifndef FORCEGROUP_H
#define FORCEGROUP_H

#include <list>
#include <vector>

#include <protomol/base/MakeableDefinition.h>
namespace ProtoMol {
  class Force;
  class MetaForce;
  class MollyForce;
  class SystemForce;
  class ExtendedForce;
  class GenericTopology;
  class Vector3DBlock;
  class ScalarStructure;
  class ReducedHessAngle;
  class ProtoMolApp;

  //________________________________________ ForceGroup

  class ForceGroup {
    // This class contains a group of forces (both system and extended).
    // It provides a convenient abstraction for integrators and other parts
    // of the system that handle user-specified groups of forces.

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ForceGroup();
    // Create an empty ForceGroup object.

    ~ForceGroup();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class ForceGroup
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    /// Evaluate all system forces in this group.
    void evaluateSystemForces(ProtoMolApp *app, Vector3DBlock *forces) const;

    /// Evaluate all extended forces in this group.
    void evaluateExtendedForces(ProtoMolApp *app, Vector3DBlock *forces) const;

    // Evaluate all MOLLY forces in this group.
    void evaluateMollyForces(GenericTopology *topo,
                             const Vector3DBlock *positions,
                             std::vector<ReducedHessAngle> *angleFilter) const;

    /// Determine if there are any forces in this group.
    bool anyForces(void) const;

    /// Determine if there are any system forces in this group.
    bool anySystemForces(void) const;

    /// Determine if there are any extended forces in this group.
    bool anyExtendedForces(void) const;

    /// Determine if there are any MOLLY forces in this group.
    bool anyMollyForces(void) const;

    /// Determine if there are any meta forces in this group.
    bool anyMetaForces(void) const;

    /// Add a new system force to this group.
    void addSystemForce(SystemForce *force);

    /// Add a new extended force to this group.
    void addExtendedForce(ExtendedForce *force);

    /// Add a new MOLLY force to this group.
    void addMollyForce(MollyForce *force);

    /// Retrieve all forces
    std::vector<Force *> getForces() const;

    /// Retrieve all forces stored in meta forces
    std::vector<Force *> getDeepMetaForces() const;

    /// Add a new meta force to this group.
    void addMetaForce(MetaForce *force);

    /// Add a new force, just cast it and add it to the correct list if
    /// non-zero ...
    void addForce(Force *force);

    void getDefinition(std::vector<MakeableDefinition> &forces) const;

    void uncache();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    std::list<SystemForce *>   mySystemForcesList;
    std::list<ExtendedForce *> myExtendedForcesList;
    std::list<MollyForce *>    myMollyForcesList;
    std::list<MetaForce *>     myMetaForcesList;
  };
  //________________________________________ INLINES

  inline bool ForceGroup::anySystemForces(void) const {
    return !mySystemForcesList.empty();
  }


  inline bool ForceGroup::anyExtendedForces(void) const {
    return !myExtendedForcesList.empty();
  }

  inline bool ForceGroup::anyMollyForces(void) const {
    return !myMollyForcesList.empty();
  }

  inline bool ForceGroup::anyMetaForces(void) const {
    return !myMetaForcesList.empty();
  }

  inline bool ForceGroup::anyForces(void) const {
    return mySystemForcesList.size() + myExtendedForcesList.size() +
           myMollyForcesList.size() + myMetaForcesList.size() > 0;
  }
}
#endif /* FORCEGROUP_H */
