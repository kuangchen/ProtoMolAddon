/* -*- c++ -*- */
#ifndef HARMDIHEDRALSYSTEMFORCE_H
#define HARMDIHEDRALSYSTEMFORCE_H

#include <protomol/force/bonded/MTorsionSystemForce.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/parallel/Parallel.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/Report.h>

#include <string>

using namespace ProtoMol::Report;
namespace ProtoMol {
  //____ HarmDihedralSystemForce

  template<class TBoundaryConditions>
  class HarmDihedralSystemForce :
    public MTorsionSystemForce<TBoundaryConditions> {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    // Default constructor
    HarmDihedralSystemForce() : k(1), myDihedral(0), myDihedralReference(0.0),
      computeOthers(false) {}

    // Constructor with parameters
    HarmDihedralSystemForce(Real kbias, int dihedral, Real dihedralReference,
                            bool other) : k(kbias), myDihedral(dihedral),
      myDihedralReference(dihedralReference), computeOthers(other) {}

    virtual ~HarmDihedralSystemForce() {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class HarmDihedralSystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    void harmCalcTorsion(const GenericTopology *topo,
                         const TBoundaryConditions &boundary,
                         const Torsion &currentTorsion,
                         const Vector3DBlock *positions, Vector3DBlock *forces,
                         Real &energy,
                         ScalarStructure *energies);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class SystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    // Avoid compiler warning/error
    using MTorsionSystemForce<TBoundaryConditions>::evaluate;
    virtual void evaluate(const GenericTopology *topo,
                          const Vector3DBlock *positions, Vector3DBlock *forces,
                          ScalarStructure *energies);

    virtual void parallelEvaluate(const GenericTopology *topo,
                                  const Vector3DBlock *positions,
                                  Vector3DBlock *forces,
                                  ScalarStructure *energies);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Force
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getKeyword() const {return "HarmDihedral";}

    virtual unsigned int numberOfBlocks(const GenericTopology *topo,
                                        const Vector3DBlock *pos);

  private:
    virtual Force *doMake(const std::vector<Value> &values) const {
      return new HarmDihedralSystemForce(values[0], values[1], values[2],
                                         values[3]);
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const {return getKeyword();}
    virtual void getParameters(std::vector<Parameter> &parameters) const {
      parameters.push_back
        (Parameter("-kbias", Value(k, ConstraintValueType::NotNegative()),
                   Text("potential bias constant")));
      parameters.push_back
        (Parameter("-dihedral",
                   Value(myDihedral, ConstraintValueType::NotNegative())));
      parameters.push_back(Parameter("-angle", Value(rtod(myDihedralReference)),
                                     Text("reference angle -180 to 180")));
      parameters.push_back
        (Parameter("-others",
                   Value(computeOthers, ConstraintValueType::NoConstraints())));
    }

  public:
    virtual void doSetParameters(std::vector<Value> values) {
      k = values[0];
      myDihedral = values[1];
      myDihedralReference = values[2];
      computeOthers = values[3];
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    // The harmonic potential biasing constant
    Real k;
    // The dihedral reference angle ID
    int myDihedral;
    // The dihedral reference angle value
    Real myDihedralReference;
    bool computeOthers;
  };

  //____ INLINES

  template<class TBoundaryConditions>
  inline void HarmDihedralSystemForce<TBoundaryConditions>::evaluate(
    const GenericTopology *topo, const Vector3DBlock *positions,
    Vector3DBlock *forces, ScalarStructure *energies) {
    const TBoundaryConditions &boundary =
      ((SemiGenericTopology<TBoundaryConditions> &)(*topo)).
        boundaryConditions;

    // Examine each dihedral
    for (unsigned int i = 0; i < topo->dihedrals.size(); i++) {
      if (computeOthers)
        calcTorsion(boundary, topo->dihedrals[i], positions, forces,
                    (*energies)[ScalarStructure::DIHEDRAL],
                    energies);
      if (static_cast<int>(i) == myDihedral)
        harmCalcTorsion(topo, boundary, topo->dihedrals[i], positions, forces,
                        (*energies)[ScalarStructure::DIHEDRAL], energies);
    }
  }

  template<class TBoundaryConditions>
  inline void HarmDihedralSystemForce<TBoundaryConditions>::parallelEvaluate(
    const GenericTopology *topo, const Vector3DBlock *positions,
    Vector3DBlock *forces, ScalarStructure *energies) {
    const TBoundaryConditions &boundary =
      (dynamic_cast<const SemiGenericTopology<TBoundaryConditions> &>(*topo)).
        boundaryConditions;

    unsigned int n = topo->dihedrals.size();
    unsigned int count = numberOfBlocks(topo, positions);

    for (unsigned int i = 0; i < count; i++)
      if (Parallel::next()) {
        int to = (n * (i + 1)) / count;
        if (to > static_cast<int>(n))
          to = n;
        int from = (n * i) / count;
        for (int j = from; j < to; j++)
          if (j != myDihedral)
            calcTorsion(boundary, topo->dihedrals[j], positions, forces,
                        (*energies)[ScalarStructure::DIHEDRAL], energies);
          else
            harmCalcTorsion(topo, boundary, topo->dihedrals[j], positions,
                            forces, (*energies)[ScalarStructure::DIHEDRAL],
                            energies);
      }
  }

  template<class TBoundaryConditions>
  inline void HarmDihedralSystemForce<TBoundaryConditions>::harmCalcTorsion(
    const GenericTopology *topo, const TBoundaryConditions &boundary,
    const Torsion & /*currTorsion*/, const Vector3DBlock *positions,
    Vector3DBlock *forces, Real &energy, ScalarStructure *energies) {
    // Calculate the energy and forces due to the harmonic potential
    // -------------------------------------------------------------
    // The dihedral angle
    Real dihedralAngle = computePhiDihedral(topo, positions, myDihedral);

    // -------------------------------------------------------------
    // The dihedral potential
    //Real k = 3.0; // kcal / (mol*K)

    Real diff = (dihedralAngle - myDihedralReference);
    if (diff < -M_PI)
      diff += 2 * M_PI;
    else if (diff > M_PI)
      diff -= 2 * M_PI;
    Real V = k * diff * diff;

    // -------------------------------------------------------------
    // The atoms in the dihedral
    int ai = topo->dihedrals[myDihedral].atom1;
    int aj = topo->dihedrals[myDihedral].atom2;
    int ak = topo->dihedrals[myDihedral].atom3;
    int al = topo->dihedrals[myDihedral].atom4;

    // -------------------------------------------------------------
    // The vector coordinates of the atoms
    Vector3D ri((*positions)[ai]);
    Vector3D rj((*positions)[aj]);
    Vector3D rk((*positions)[ak]);
    Vector3D rl((*positions)[al]);
    // -------------------------------------------------------------
    // Vectors between atoms (rij = ri - rj)
    Vector3D rij = boundary.minimalDifference(rj, ri);
    Vector3D rkj = boundary.minimalDifference(rj, rk);
    Vector3D rkl = boundary.minimalDifference(rl, rk);
    // -------------------------------------------------------------
    // Normals
    Vector3D m = rij.cross(rkj);
    Vector3D n = rkj.cross(rkl);
    // -------------------------------------------------------------
    // The derivate of V with respect to the dihedral angle
    Real dVdPhi = 2 * (k) * diff;

    // -------------------------------------------------------------
    // Miscellaneous quantities needed to compute the forces
    Real rkj_norm = rkj.norm();
    Real rkj_normsq = rkj.normSquared();
    Real m_normsq = m.normSquared();
    Real n_normsq = n.normSquared();
    Real rij_dot_rkj = rij.dot(rkj);
    Real rkl_dot_rkj = rkl.dot(rkj);
    // -------------------------------------------------------------
    // Forces on the atoms
    // Atom i
    Vector3D fi = m * (-dVdPhi * rkj_norm / m_normsq);
    // Atom l
    Vector3D fl = n * (dVdPhi * rkj_norm / n_normsq);
    // Atom j
    Vector3D fj = fi *
                  (-1 + rij_dot_rkj /
                   rkj_normsq) - fl * (rkl_dot_rkj / rkj_normsq);
    // Atom k
    Vector3D fk = -(fi + fj + fl);

    // *************************************************************
    // -------------------------------------------------------------
    // Update the dihedral system energy
    energy += V;
    // Update the atom forces
    (*forces)[ai] += fi;
    (*forces)[aj] += fj;
    (*forces)[ak] += fk;
    (*forces)[al] += fl;
    // -------------------------------------------------------------
    // Compute the virial energy
    if (energies->virial()) {
      Vector3D f1 = fi;
      Vector3D f2 = fi + fj;
      Vector3D f3 = -fl;
      Vector3D r12 = rij;
      Vector3D r23 = -rkj;
      Vector3D r34 = rkl;

      Real xy = f1.c[0] * r12.c[1] + f2.c[0] * r23.c[1] + f3.c[0] * r34.c[1];
      Real xz = f1.c[0] * r12.c[2] + f2.c[0] * r23.c[2] + f3.c[0] * r34.c[2];
      Real yz = f1.c[1] * r12.c[2] + f2.c[1] * r23.c[2] + f3.c[1] * r34.c[2];

      (*energies)[ScalarStructure::VIRIALXX] += f1.c[0] * r12.c[0] + f2.c[0] * r23.c[0] +
                                                f3.c[0] * r34.c[0];
      (*energies)[ScalarStructure::VIRIALXY] += xy;
      (*energies)[ScalarStructure::VIRIALXZ] += xz;
      (*energies)[ScalarStructure::VIRIALYX] += xy;
      (*energies)[ScalarStructure::VIRIALYY] += f1.c[1] * r12.c[1] + f2.c[1] * r23.c[1] +
                                                f3.c[1] * r34.c[1];
      (*energies)[ScalarStructure::VIRIALYZ] += yz;
      (*energies)[ScalarStructure::VIRIALZX] += xz;
      (*energies)[ScalarStructure::VIRIALZY] += yz;
      (*energies)[ScalarStructure::VIRIALZZ] += f1.c[2] * r12.c[2] + f2.c[2] * r23.c[2] +
                                                f3.c[2] * r34.c[2];
    }
  }

  template<class TBoundaryConditions>
  inline unsigned int HarmDihedralSystemForce<TBoundaryConditions>::
    numberOfBlocks(const GenericTopology *topo, const Vector3DBlock *) {
    return std::min(Parallel::getAvailableNum(),
                    static_cast<int>(topo->dihedrals.size()));
  }
}
#endif /* HARMDIHEDRALSYSTEMFORCE_H */
