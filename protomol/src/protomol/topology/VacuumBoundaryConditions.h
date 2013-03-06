/* -*- c++ -*- */
#ifndef VACUUMBOUNDARYCONDITIONS_H
#define VACUUMBOUNDARYCONDITIONS_H

#include <protomol/config/Parameter.h>
#include <protomol/type/Vector3DBlock.h>
#include <string>
#include <iostream>
using namespace std;

namespace ProtoMol {
  //________________________________________ VacuumBoundaryConditions

  /**
   * These are the vacuum (non-periodic or normal/standard) boundary conditions.
   */
  class VacuumBoundaryConditions {
  public:

    /// Set method for the dimensions of the (original) simulation box.
    void set(const Vector3D &, const Vector3D &, const Vector3D &,
             const Vector3D &) const {};

    /// Perform a minimal-image subtraction.
    Vector3D minimalDifference(const Vector3D &c1, const Vector3D &c2) const {
      Vector3D diff(c2);
      diff -= c1;
      return diff;
    }

    // Perform a minimal-image subtraction and computes the squared distance
    Vector3D minimalDifference(const Vector3D &c1, const Vector3D &c2,
                               Real &distSquared) const {
      Vector3D diff(c2);
      diff -= c1;
      distSquared = diff.normSquared();
      return diff;
    }
    /// Find the position in the basis/original cell/image.
    Vector3D minimalPosition(const Vector3D &c) const {return c;}
    /// Find the lattice vector difference between two positions
    Vector3D minimalTranslationDifference(const Vector3D &,
                                          const Vector3D &) const {
      return Vector3D(0.0, 0.0, 0.0);
    }
    /// Find the lattice translation relative to the original cell/image
    Vector3D minimalTranslationPosition(const Vector3D &) const {
      return Vector3D(0.0, 0.0, 0.0);
    }


    /// basis/unit vector e1
    Vector3D e1()        const;
    /// basis/unit vector e2
    Vector3D e2()        const;
    /// basis/unit vector e3
    Vector3D e3()        const;
    /// inverse basis/unit vector e1
    Vector3D e1r()       const;
    /// inverse basis/unit vector e2
    Vector3D e2r()       const;
    /// inverse basis/unit vector e3
    Vector3D e3r()       const;
    /// origin of the minimal image
    Vector3D origin()    const;
    /// minimal corner of the bounding box of the minimal image/cell
    Vector3D getMin()    const;
    /// maximal corner of the bounding box of the minimal image/cell
    Vector3D getMax()    const;
    Real getVolume()     const;
    bool isOrthogonal()  const {return true;};

    /// Boolean's defining boundary conditions
    enum {PERIODIC = 0, VACUUM = 1};

    /**
     * Builds the lattice vector for a given cutoff. Since vacuum boundary 
     * conditions have no replications and (0,0,0) always is excluded, there 
     * is no need for any lattice vectors.
     */
    std::vector<Vector3D> buildLatticeVectors(Real cutoff) const;

    const std::string &getKeyword() const {return keyword;}
    void getParameters(std::vector<Parameter> &) const {};
    static VacuumBoundaryConditions make(std::vector<Value> values);
    static unsigned int getParameterSize() {return 0;}
    std::vector<Parameter> getDefaults(const Vector3DBlock &) const {
      return std::vector<Parameter>();
    }
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
  private:
  };

  //________________________________________ INLINES
}
#endif /* VACUUMBOUNDARYCONDITIONS_H */
