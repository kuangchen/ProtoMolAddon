/*  -*- c++ -*-  */
#ifndef ONEATOMCONSTRAINTS_H
#define ONEATOMCONSTRAINTS_H

#include <protomol/topology/GenericTopology.h>

namespace ProtoMol {
  /**
   * Retrieves the cutoff^2 from either the switching function or the
   * potential. TNonbondedForce::CUTOFF should be true
   * if a look up table is used.
   */
  template<bool flag> struct Cutoff {
    template<class TS, class TF>
    static Real cutoff(const TS &swf, const TF &) {return swf.cutoffSquared();}
  };
  template<> struct Cutoff<true> {
    template<class TS, class TF>
    static Real cutoff(const TS &, const TF &nf) {return nf.cutoffSquared();}
  };

  /*
     Constraints for OneAtomPair's templates and 1-body
     forces/potentials. The constraint is checked at beginning of
     doOneAtomPair(), which makes pretty efficient. If check() fails
     neither the distance,    exclusions, potential nor switching
     function are computed. getPrefixId() and getPostfixId() are added
     before and after the potential keyword in order to distinguish
     between a potential with and one without constraint.
   */

  /**
     Default constraint of OneAtomPair's templates, for all atoms pairs true
   */
  struct NoConstraint {
    enum {PRE_CHECK = 0, POST_CHECK = 0};
    static bool check(const GenericTopology *, int, int) {return true;}

    static bool check(const GenericTopology *, int) {return true;}

    static std::string getPrefixId() {return "";}

    static std::string getPostfixId() {return "";}

    static void check(const GenericTopology *, int, int, const Vector3D &, Real,
                      const Vector3D &) {}
  };

  /**
     HBond constraint of OneAtomPair's templates, only pairs with different
     mass are considered
   */
  struct HBondConstraint {
    enum {PRE_CHECK = 1, POST_CHECK = 0};
    static bool check(const GenericTopology *topo, int i, int j) {
      return !equal(topo->atoms[i].scaledMass, topo->atoms[j].scaledMass);
    }

    static std::string getPrefixId() {return "HBond";}

    static std::string getPostfixId() {return "";}

    static void check(const GenericTopology *, int, int, const Vector3D &, Real,
                      const Vector3D &) {}
  };

  /**
     Constraint of OneAtomPair's templates, only pairs on the same molecule
     are considered
   */
  struct SameMoleculeConstraint {
    enum {PRE_CHECK = 1, POST_CHECK = 0};
    static bool check(const GenericTopology *topo, int i,
                      int j) {return topo->atoms[i].molecule ==
                                     topo->atoms[j].molecule;}

    static std::string getPrefixId() {return "";}

    static std::string getPostfixId() {return "SameMol";}

    static void check(const GenericTopology *, int, int, const Vector3D &, Real,
                      const Vector3D &) {}
  };

  /**
     Constraint of OneAtomPair's templates, only pairs off different moleculea
     are considered
   */
  struct NotSameMoleculeConstraint {
    enum {PRE_CHECK = 1, POST_CHECK = 0};
    static bool check(const GenericTopology *topo, int i,
                      int j) {return topo->atoms[i].molecule !=
                                     topo->atoms[j].molecule;}

    static std::string getPrefixId() {return "";}

    static std::string getPostfixId() {return "NotSameMol";}

    static void check(const GenericTopology *, int, int, const Vector3D &, Real,
                      const Vector3D &) {}
  };

  /**
     Constraint of OneAtomPair's templates, print pair
   */
  struct DebugPreConstraint {
    enum {PRE_CHECK = 1, POST_CHECK = 0};
    static bool check(const GenericTopology *, int i, int j) {
      Report::report << "P(" <<
      std::min(i, j) << "," << std::max(i, j) << ")" << Report::endr;
      return true;
    }

    static bool check(const GenericTopology *, int) {return true;}

    static std::string getPrefixId() {return "";}

    static std::string getPostfixId() {return "DebugPre";}

    static void check(const GenericTopology *, int, int, const Vector3D &, Real,
                      const Vector3D &) {}
  };

  struct DebugPostConstraint {
    enum {PRE_CHECK = 0, POST_CHECK = 1};
    static bool check(const GenericTopology *, int, int) {return true;}

    static bool check(const GenericTopology *, int) {return true;}

    static std::string getPrefixId() {return "";}

    static std::string getPostfixId() {return "DebugPost";}

    static void check(const GenericTopology *, int i, int j,
                      const Vector3D &diff, Real energy,
                      const Vector3D &force) {
      Report::report.precision(15);
      Report::report << "P(" <<
      std::min(i,
               j) << "," <<
      std::max(i,
               j) << ") " <<
      (i <
       j ? diff : -diff) << ", " << energy << ", " <<
      (i < j ? force : -force) << Report::endr;
    }
  };
}
#endif /* ONEATOMPAIRCONSTRAINTS_H */
