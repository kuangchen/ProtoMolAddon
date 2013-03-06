/* -*- c++ -*- */
#ifndef PERODICBOUNDARYCONDITIONS_H
#define PERODICBOUNDARYCONDITIONS_H

#include <string>

#include <protomol/config/Parameter.h>
#include <protomol/type/Vector3DBlock.h>

namespace ProtoMol {
  //________________________________________ PeriodicBoundaryConditions

  /**
   * Implements periodic boundary conditions, defining how we measure distances
   * and accounting the wrapping-around effect.
   * The class use a couple of shorts cut's to avoid rint and to many div's and
   * mul's.
   */
  class PeriodicBoundaryConditions {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    PeriodicBoundaryConditions();
    PeriodicBoundaryConditions(const Vector3D &e1, const Vector3D &e2,
                               const Vector3D &e3, const Vector3D &origin);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class PeriodicBoundaryConditions
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    /// Set method for the dimensions of the (original) simulation box.
    void set(const Vector3D &e1, const Vector3D &e2, const Vector3D &e3,
             const Vector3D &origin);

    /// Perform a minimal-image subtraction.
    Vector3D minimalDifference(const Vector3D &c1, const Vector3D &c2) const {
      Vector3D diff(c2);
      diff -= c1;
      if (myOrthogonal) {
        if (diff.normSquared() > myD) {
          // diff not small engouh, we do have to do more ...
          // ... may be it's just a single wrapping
          if (diff.c[0] < -myH2.c[0])
            diff.c[0] += myH.c[0];
          else if (diff.c[0] > myH2.c[0])
            diff.c[0] -= myH.c[0];
          if (diff.c[1] < -myH2.c[1])
            diff.c[1] += myH.c[1];
          else if (diff.c[1] > myH2.c[1])
            diff.c[1] -= myH.c[1];
          if (diff.c[2] < -myH2.c[2])
            diff.c[2] += myH.c[2];
          else if (diff.c[2] > myH2.c[2])
            diff.c[2] -= myH.c[2];
          if (diff.normSquared() > myD) {
            // ... distance was pretty big, hu? ...
            diff.c[0] -= myE1.c[0] * rint(myE1r.c[0] * diff.c[0]);
            diff.c[1] -= myE2.c[1] * rint(myE2r.c[1] * diff.c[1]);
            diff.c[2] -= myE3.c[2] * rint(myE3r.c[2] * diff.c[2]);
          }
        }
      } else
        // ... really difficult!
        diff -=
          Vector3D(myE1 * rint(myE1r.dot(diff)) + myE2 * rint(myE2r.dot(
                diff)) + myE3 * rint(myE3r.dot(diff)));
      return diff;
    }

    // Perform a minimal-image subtraction and computes the squared distance
    Vector3D minimalDifference(const Vector3D &c1, const Vector3D &c2,
                               Real &distSquared) const {
      Vector3D diff(c2);
      diff -= c1;
      if (myOrthogonal) {
        distSquared = diff.normSquared();
        if (distSquared > myD) {
          // diff not small engouh, we do have to do more ...
          // ... may be it's just a single wrapping
          if (diff.c[0] < -myH2.c[0])
            diff.c[0] += myH.c[0];
          else if (diff.c[0] > myH2.c[0])
            diff.c[0] -= myH.c[0];
          if (diff.c[1] < -myH2.c[1])
            diff.c[1] += myH.c[1];
          else if (diff.c[1] > myH2.c[1])
            diff.c[1] -= myH.c[1];
          if (diff.c[2] < -myH2.c[2])
            diff.c[2] += myH.c[2];
          else if (diff.c[2] > myH2.c[2])
            diff.c[2] -= myH.c[2];
          distSquared = diff.normSquared();
          if (distSquared > myD) {
            // ... distance was pretty big, hu? ...
            diff.c[0] -= myE1.c[0] * rint(myE1r.c[0] * diff.c[0]);
            diff.c[1] -= myE2.c[1] * rint(myE2r.c[1] * diff.c[1]);
            diff.c[2] -= myE3.c[2] * rint(myE3r.c[2] * diff.c[2]);
            distSquared = diff.normSquared();
          }
        }
      } else {
        // ... really difficult!
        diff -=
          Vector3D(myE1 * rint(myE1r.dot(diff)) + myE2 * rint(myE2r.dot(
                diff)) + myE3 * rint(myE3r.dot(diff)));
        distSquared = diff.normSquared();
      }
      return diff;
    }

    /// Find the position in the basis/original cell/image.
    Vector3D minimalPosition(const Vector3D &c) const {
      Vector3D diff(c);
      diff -= myOrigin;
      if (myOrthogonal)
        if (diff.normSquared() <= myD)
          // diff so small, we do not have to do more ...
          return diff;
        else {
          // ... may be it's just a single wrapping
          if (diff.c[0] < -myH2.c[0])
            diff.c[0] += myH.c[0];
          else if (diff.c[0] > myH2.c[0])
            diff.c[0] -= myH.c[0];
          if (diff.c[1] < -myH2.c[1])
            diff.c[1] += myH.c[1];
          else if (diff.c[1] > myH2.c[1])
            diff.c[1] -= myH.c[1];
          if (diff.c[2] < -myH2.c[2])
            diff.c[2] += myH.c[2];
          else if (diff.c[2] > myH2.c[2])
            diff.c[2] -= myH.c[2];
          if (diff.normSquared() <= myD)
            return diff;
          else
            // ... distance was pretty big, hu? ...
            return Vector3D(diff - Vector3D(myE1.c[0] * rint(myE1r.c[0] * diff.c[0]),
                myE2.c[1] * rint(myE2r.c[1] * diff.c[1]),
                myE3.c[2] * rint(myE3r.c[2] * diff.c[2])));
        }
      else
        // ... really difficult!
        return Vector3D(diff - myE1 * rint(myE1r.dot(diff))
          - myE2 * rint(myE2r.dot(diff))
          - myE3 * rint(myE3r.dot(diff)));
    }

    /// Find the lattice vector difference between two positions
    Vector3D minimalTranslationDifference(const Vector3D &c1,
                                          const Vector3D &c2) const {
      Vector3D diff(c2 - c1);
      if (myOrthogonal)
        if (diff.normSquared() <= myD)
          // diff so small, we do not have to do more ...
          return Vector3D(0.0, 0.0, 0.0);
        else
          return Vector3D(myE1.c[0] * rint(myE1r.c[0] * diff.c[0]),
            myE2.c[1] * rint(myE2r.c[1] * diff.c[1]),
            myE3.c[2] * rint(myE3r.c[2] * diff.c[2]));
      else
        // ... really difficult!
        return Vector3D(myE1 * rint(myE1r.dot(diff))
          + myE2 * rint(myE2r.dot(diff))
          + myE3 * rint(myE3r.dot(diff)));
    }

    /// Find the lattice translation relative to the original cell/image
    Vector3D minimalTranslationPosition(const Vector3D &c) const {
      Vector3D diff(c - myOrigin);
      if (myOrthogonal)
        if (diff.normSquared() <= myD)
          // diff so small, we do not have to do more ...
          return Vector3D(0.0, 0.0, 0.0);
        else
          return Vector3D(myE1.c[0] * rint(myE1r.c[0] * diff.c[0]),
            myE2.c[1] * rint(myE2r.c[1] * diff.c[1]),
            myE3.c[2] * rint(myE3r.c[2] * diff.c[2]));
      else
        // ... really difficult!
        return Vector3D(myE1 * rint(myE1r.dot(diff))
          + myE2 * rint(myE2r.dot(diff))
          + myE3 * rint(myE3r.dot(diff)));

    }

    /// basis/unit vector e1
    const Vector3D &e1()     const {return myE1;}
    /// basis/unit vector e2
    const Vector3D &e2()     const {return myE2;}
    /// basis/unit vector e3
    const Vector3D &e3()     const {return myE3;}
    /// inverse basis/unit vector e1
    const Vector3D &e1r()    const {return myE1r;}
    /// inverse basis/unit vector e2
    const Vector3D &e2r()    const {return myE2r;}
    /// inverse basis/unit vector e3
    const Vector3D &e3r()    const {return myE3r;}
    /// origin of the minimal image
    const Vector3D &origin() const {return myOrigin;}
    /// minimal corner of the bounding box of the minimal image/cell
    const Vector3D &getMin() const {return myMin;}
    /// maximal corner of the bounding box of the minimal image/cell
    const Vector3D &getMax() const {return myMax;}
    Real getVolume()         const {return myV;};
    bool isOrthogonal()      const {return myOrthogonal;};

    /// Boolean's defining boundary conditions
    enum {PERIODIC = 1, VACUUM = 0};

    /// Builds the lattice vector for a given cutoff. (0,0,0) is not included.
    std::vector<Vector3D> buildLatticeVectors(Real cutoff) const;

    /// Returns the keyword of the boundary conditions
    const std::string &getKeyword() const {return keyword;}
    void getParameters(std::vector<Parameter> &parameters) const;
    static PeriodicBoundaryConditions make(std::vector<Value> values);
    static unsigned int getParameterSize() {return 4;}

    /// Returns possible default values for the parameters based on the
    /// positions
    std::vector<Parameter> getDefaults(const Vector3DBlock &positions) const;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
  private:
    Vector3D myE1;
    Vector3D myE2;
    Vector3D myE3;
    Vector3D myE1r;
    Vector3D myE2r;
    Vector3D myE3r;
    Vector3D myOrigin;
    Vector3D myMin;
    Vector3D myMax;

    Real myDX;
    Real myDY;
    Real myDZ;
    /// maximal distance between two positions where plain subtraction if safe
    Real myD;  
    Vector3D myH;
    Vector3D myH2;

    Real myV;
    bool myOrthogonal;
  };
}
#endif /* PERODICBOUNDARYCONDITIONS_H */
