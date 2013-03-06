/* -*- c++ -*- */
#ifndef SEMIGENERICTOPOLOGY_H
#define SEMIGENERICTOPOLOGY_H

#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>

namespace ProtoMol {
  //________________________________________ SemiGenericTopology
  /**
   * Introducing a semi-generic topology for the forces which do not
   * depend on the cell manager
   */
  template<class TBoundaryConditions>
  class SemiGenericTopology : public GenericTopology {
  public:
    typedef TBoundaryConditions BoundaryConditions;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    SemiGenericTopology() : GenericTopology() {}
    SemiGenericTopology(Real csf, const ExclusionType &e,
                        const TBoundaryConditions &b) :
      GenericTopology(csf, e), boundaryConditions(b) {}
    virtual ~SemiGenericTopology() {};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class GenericTopology
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    /// Returns possible default values for the parameters based on the
    /// positions
    virtual std::vector<Parameter> getDefaults(const Vector3DBlock &positions)
    const {return boundaryConditions.getDefaults(positions);}
    /// Perform a minimal-image subtraction.
    virtual Vector3D minimalDifference(const Vector3D &c1,
                                       const Vector3D &c2) const {
      return boundaryConditions.minimalDifference(c1, c2);
    }
    /// Find the position in the basis/original cell/image.
    virtual Vector3D minimalPosition(const Vector3D &c) const {
      return Vector3D(boundaryConditions.minimalPosition(c));
    }
    virtual Real getVolume() const {return boundaryConditions.getVolume();}

    virtual Real getVolume(const Vector3DBlock &positions) const {
      Real v = boundaryConditions.getVolume();
      if (v >= Constant::MAXREAL) {
        Vector3D a, b;
        positions.boundingbox(a, b);
        v = (b.c[0] - a.c[0]) * (b.c[1] - a.c[1]) * (b.c[2] - a.c[2]);
      }
      return v;
    }

    /// rescales the volume of the system cell/minimal image
    virtual void rescaleVolume(Real fac) {
      Real f = pow(fac, 1.0 / 3.0);
      boundaryConditions.set(boundaryConditions.e1() * f,
        boundaryConditions.e2() * f,
        boundaryConditions.e3() * f,
        boundaryConditions.origin());
    }

    /// computes the minimal bounding box of the positions taking into account
    /// of the boundary conditions
    virtual void getBoundingbox(const Vector3DBlock &positions,
                                Vector3D &minbb, Vector3D &maxbb) const {
      if (positions.size() <= 0) {
        minbb =
          Vector3D(Constant::MAXREAL, Constant::MAXREAL, Constant::MAXREAL);
        maxbb =
          Vector3D(-Constant::MAXREAL, -Constant::MAXREAL, -Constant::MAXREAL);
        return;
      }

      Vector3D origin(boundaryConditions.origin());

      minbb = boundaryConditions.minimalPosition(positions[0]) + origin;
      maxbb = minbb;

      const unsigned int count = positions.size();
      for (unsigned int i = 1; i < count; ++i) {
        Vector3D pos(boundaryConditions.minimalPosition(positions[i]) + origin);
        if (pos.c[0] < minbb.c[0]) minbb.c[0] = pos.c[0];
        else if (pos.c[0] > maxbb.c[0]) maxbb.c[0] = pos.c[0];
        if (pos.c[1] < minbb.c[1]) minbb.c[1] = pos.c[1];
        else if (pos.c[1] > maxbb.c[1]) maxbb.c[1] = pos.c[1];
        if (pos.c[2] < minbb.c[2]) minbb.c[2] = pos.c[2];
        else if (pos.c[2] > maxbb.c[2]) maxbb.c[2] = pos.c[2];
      }
    }

    /// returns the bounding box of the boundary conditions
    virtual void getBoundaryConditionsBox(Vector3D &minbb,
                                          Vector3D &maxbb) const {
      minbb = boundaryConditions.getMin();
      maxbb = boundaryConditions.getMax();
    }

    /**
     * checks whether the plain distances between all atoms on each molecule are
     * the same as with applying the boundary conditions. If true, there is no 
     * molecule, which is wrapped around.
     */
    virtual bool checkMoleculePairDistances(const Vector3DBlock &positions)
    const {
      for (unsigned int i = 0; i < this->molecules.size(); ++i)
        for (unsigned int j = 0; j < this->molecules[i].pairs.size(); ++j) {
          Vector3D c1(positions[this->molecules[i].pairs[j].first]);
          Vector3D c2(positions[this->molecules[i].pairs[j].second]);
          if (boundaryConditions.minimalTranslationDifference
              (c1, c2).normSquared() > Constant::EPSILON)
            return false;
        }

      return true;
    }

    /**
     * translates all positions of a each molecule such that the plain 
     * distances between them yields and that the center of mass is inside the 
     * minimal image/basis cell. It may not work for molecules which are bigger
     * then half of the simulation cell.
     */
    virtual void minimalImage(Vector3DBlock &positions) {
      for (unsigned int i = 0; i < this->molecules.size(); ++i) {
        for (unsigned int j = 0; j < this->molecules[i].pairs.size(); ++j) {
          Vector3D c1(positions[this->molecules[i].pairs[j].first]);
          Vector3D c2(positions[this->molecules[i].pairs[j].second]);
          Vector3D d(boundaryConditions.minimalTranslationDifference(c1, c2));

          if (d.normSquared() > Constant::EPSILON)
            positions[this->molecules[i].pairs[j].second] -= d;
        }

        Vector3D center(molecularCenterOfMass(this->molecules[i].atoms,
                                              &positions, this));
        Vector3D d(boundaryConditions.minimalTranslationPosition(center));

        if (d.normSquared() > Constant::EPSILON) {
          center -= d;
          for (unsigned int j = 0; j < this->molecules[i].atoms.size(); ++j)
            positions[this->molecules[i].atoms[j]] -= d;
        }

        this->molecules[i].position = center;
      }
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    BoundaryConditions boundaryConditions;
  };
}
#endif /* SEMIGENERICTOPOLOGY_H */
