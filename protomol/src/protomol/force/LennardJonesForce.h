/* -*- c++ -*- */
#ifndef LENNARDJONESFORCE_H
#define LENNARDJONESFORCE_H

#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/LennardJonesParameters.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/topology/ExclusionTable.h>
#include <protomol/config/Parameter.h>
#include <string>
#include <protomol/base/Report.h>
using namespace ProtoMol::Report;

namespace ProtoMol {
  //____ LennardJonesForce
  class LennardJonesForce {
  public:
    enum {DIST_R2 = 1};
    enum {CUTOFF = 0};
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class LennardJonesForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void operator()(Real &energy, Real &force,
                    Real /*distSquared*/, Real rDistSquared, const Vector3D &,
                    const GenericTopology *topo, int atom1, int atom2,
                    ExclusionClass excl) const {
#ifdef DEBUG_LJ_DISTANCE_CHECK
      if (rDistSquared < 1.0 / 0.25)
        report << warning << "LennardJonesForce::operator(): atom " << atom1 <<
        " and " << atom2 << " get closer"
                << " than 0.5 angstroms (=" <<
        sqrt(1.0 / rDistSquared) << ")." << endr;
#endif
      const LennardJonesParameters &
      params = topo->lennardJonesParameters(topo->atoms[atom1].type,
                                            topo->atoms[atom2].type);

      Real A, B;
      if (excl != EXCLUSION_MODIFIED) {
        A = params.A;
        B = params.B;
      } else {
        A = params.A14;
        B = params.B14;
      }

      // Fast LJ computation
      //Real r2   = rDistSquared;
      Real r6 = rDistSquared * rDistSquared * rDistSquared;
      Real r12 = r6 * r6;
      Real r6B = B * r6;
      Real r12A = A * r12;
      energy = r12A - r6B;
      force = 12.0 * r12A * rDistSquared - 6.0 * r6B * rDistSquared;
    }

    static void accumulateEnergy(ScalarStructure *energies, Real energy) {
      (*energies)[ScalarStructure::LENNARDJONES] += energy;
    }

    static Real getEnergy(const ScalarStructure *energies) {
      return (*energies)[ScalarStructure::LENNARDJONES];
    }

    // Parsing
    static std::string getId() {return keyword;}
    static unsigned int getParameterSize() {return 0;}
    void getParameters(std::vector<Parameter> &) const {}

    static LennardJonesForce make(const std::vector<Value> &) {
      return LennardJonesForce();
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
  private:
  };
  //____ INLINES
}

#endif /* LENNARDJONESFORCE_H */
