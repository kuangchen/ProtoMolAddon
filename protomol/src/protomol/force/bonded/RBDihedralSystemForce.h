/* -*- c++ -*- */
#ifndef RBDIHEDRALSYSTEMFORCE_H
#define RBDIHEDRALSYSTEMFORCE_H

#include <protomol/force/system/SystemForce.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/parallel/Parallel.h>
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/topology/RBTorsion.h>

#include <string>

namespace ProtoMol {
  //____ RBDihedralSystemForce

  template<class TBoundaryConditions>
  class RBDihedralSystemForce : public SystemForce {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual ~RBDihedralSystemForce() {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class SystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    using SystemForce::evaluate; // Avoid compiler warning/error
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
    virtual std::string getKeyword() const {return "RBDihedral";}

    virtual unsigned int numberOfBlocks(const GenericTopology *topo,
                                        const Vector3DBlock *pos);

  private:
    virtual Force *doMake(const std::vector<Value> &) const {
      return new RBDihedralSystemForce();
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const {return getKeyword();}
    virtual void getParameters(std::vector<Parameter> &) const {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My methods
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    void calcRBTorsion(const TBoundaryConditions &boundary, const RBTorsion &currRBTorsion,
              const Vector3DBlock *positions, Vector3DBlock *forces,
              Real &energy, ScalarStructure *energies ) {

        int a1 = currRBTorsion.atom1;
        int a2 = currRBTorsion.atom2;
        int a3 = currRBTorsion.atom3;
        int a4 = currRBTorsion.atom4;

        Vector3D r12(boundary.minimalDifference((*positions)[a2],
                                                (*positions)[a1]));
        Vector3D r23(boundary.minimalDifference((*positions)[a3],
                                                (*positions)[a2]));
        Vector3D r34(boundary.minimalDifference((*positions)[a4],
                                                (*positions)[a3]));

        // Cross product of r12 and r23, represents the plane shared by these two 
        // vectors
        Vector3D a(r12.cross(r23));
        // Cross product of r23 and r34, represents the plane shared by these two 
        // vectors
        Vector3D b(r23.cross(r34));
        // Cross product of r23 and A, represents the plane shared by these two 
        // vectors
        Vector3D c(r23.cross(a));

        // 1/length of Vector A, B and C
        Real ra = 1.0 / a.norm();
        Real rb = 1.0 / b.norm();
        Real rc = 1.0 / c.norm();

        // Normalize A,B and C
        a *= ra;
        b *= rb;
        c *= rc;

        // Calculate phi
        Real cosPhi = a.dot(b);
        Real sinPhi = c.dot(b);
        Real phi = -atan2(sinPhi, cosPhi);

        Real dpotdphi = 0.;
        Real Cn[6] = {currRBTorsion.C0, currRBTorsion.C1, currRBTorsion.C2, 
                      currRBTorsion.C3, currRBTorsion.C4, currRBTorsion.C5 };
        Real cosPsi = cos(phi - M_PI);
        Real sinPsi = sin(phi - M_PI);
        Real cosNm1 = 1.;

        energy += Cn[0];

        for(int i=1; i<6; i++) {

          dpotdphi -= (Real)i * Cn[i] * sinPsi * cosNm1;

          cosNm1 *= cosPsi;

          // Add energy
          energy += Cn[i] * cosNm1;

        }

        // To prevent potential singularities, if abs(sinPhi) <= 0.1, then
        // use another method of calculating the gradient.
        Vector3D f1, f2, f3;
        if (fabs(sinPhi) > 0.1) {
          //  use the sin version to avoid 1/cos terms

          Vector3D dcosdA((a * cosPhi - b) * ra);
          Vector3D dcosdB((b * cosPhi - a) * rb);

          Real k1 = dpotdphi / sinPhi;

          f1.c[0] = k1 * (r23.c[1] * dcosdA.c[2] - r23.c[2] * dcosdA.c[1]);
          f1.c[1] = k1 * (r23.c[2] * dcosdA.c[0] - r23.c[0] * dcosdA.c[2]);
          f1.c[2] = k1 * (r23.c[0] * dcosdA.c[1] - r23.c[1] * dcosdA.c[0]);

          f3.c[0] = k1 * (r23.c[2] * dcosdB.c[1] - r23.c[1] * dcosdB.c[2]);
          f3.c[1] = k1 * (r23.c[0] * dcosdB.c[2] - r23.c[2] * dcosdB.c[0]);
          f3.c[2] = k1 * (r23.c[1] * dcosdB.c[0] - r23.c[0] * dcosdB.c[1]);

          f2.c[0] = k1 *
                 (r12.c[2] * dcosdA.c[1] - r12.c[1] * dcosdA.c[2] + r34.c[1] * dcosdB.c[2] - r34.c[2] *
                  dcosdB.c[1]);
          f2.c[1] = k1 *
                 (r12.c[0] * dcosdA.c[2] - r12.c[2] * dcosdA.c[0] + r34.c[2] * dcosdB.c[0] - r34.c[0] *
                  dcosdB.c[2]);
          f2.c[2] = k1 *
                 (r12.c[1] * dcosdA.c[0] - r12.c[0] * dcosdA.c[1] + r34.c[0] * dcosdB.c[1] - r34.c[1] *
                  dcosdB.c[0]);
        } else {
          //  This angle is closer to 0 or 180 than it is to
          //  90, so use the cos version to avoid 1/sin terms

          Vector3D dsindC((c * sinPhi - b) * rc);
          Vector3D dsindB((b * sinPhi - c) * rb);

          Real k1 = -dpotdphi / cosPhi;

          f1.c[0] = k1 *
                 ((r23.c[1] * r23.c[1] + r23.c[2] *
                   r23.c[2]) * dsindC.c[0] - r23.c[0] * r23.c[1] * dsindC.c[1] - r23.c[0] * r23.c[2] *
                  dsindC.c[2]);
          f1.c[1] = k1 *
                 ((r23.c[2] * r23.c[2] + r23.c[0] *
                   r23.c[0]) * dsindC.c[1] - r23.c[1] * r23.c[2] * dsindC.c[2] - r23.c[1] * r23.c[0] *
                  dsindC.c[0]);
          f1.c[2] = k1 *
                 ((r23.c[0] * r23.c[0] + r23.c[1] *
                   r23.c[1]) * dsindC.c[2] - r23.c[2] * r23.c[0] * dsindC.c[0] - r23.c[2] * r23.c[1] *
                  dsindC.c[1]);

          f3 = dsindB.cross(r23) * k1;

          f2.c[0] = k1 *
                 (-(r23.c[1] * r12.c[1] + r23.c[2] *
                    r12.c[2]) * dsindC.c[0] +
                  (2.0 * r23.c[0] * r12.c[1] - r12.c[0] * r23.c[1]) * dsindC.c[1]
                  + (2.0 * r23.c[0] * r12.c[2] - r12.c[0] *
                     r23.c[2]) * dsindC.c[2] + dsindB.c[2] * r34.c[1] - dsindB.c[1] * r34.c[2]);
          f2.c[1] = k1 *
                 (-(r23.c[2] * r12.c[2] + r23.c[0] *
                    r12.c[0]) * dsindC.c[1] +
                  (2.0 * r23.c[1] * r12.c[2] - r12.c[1] * r23.c[2]) * dsindC.c[2]
                  + (2.0 * r23.c[1] * r12.c[0] - r12.c[1] *
                     r23.c[0]) * dsindC.c[0] + dsindB.c[0] * r34.c[2] - dsindB.c[2] * r34.c[0]);
          f2.c[2] = k1 *
                 (-(r23.c[0] * r12.c[0] + r23.c[1] *
                    r12.c[1]) * dsindC.c[2] +
                  (2.0 * r23.c[2] * r12.c[0] - r12.c[2] * r23.c[0]) * dsindC.c[0]
                  + (2.0 * r23.c[2] * r12.c[1] - r12.c[2] *
                     r23.c[1]) * dsindC.c[1] + dsindB.c[1] * r34.c[0] - dsindB.c[0] * r34.c[1]);
        }
        (*forces)[a1] += f1;
        (*forces)[a2] += f2 - f1;
        (*forces)[a3] += f3 - f2;
        (*forces)[a4] -= f3;

        // Add virial
        if (energies->virial()) {
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

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
  };

  //____ INLINES

  template<class TBoundaryConditions>
  inline void RBDihedralSystemForce<TBoundaryConditions>::
  evaluate(const GenericTopology *topo, const Vector3DBlock *positions,
           Vector3DBlock *forces, ScalarStructure *energies) {

    const TBoundaryConditions &boundary =
      ((SemiGenericTopology<TBoundaryConditions> &)(*topo)).
        boundaryConditions;

    unsigned int n = topo->rb_dihedrals.size();

    for (unsigned int i = 0; i < n; i++) {
      calcRBTorsion(boundary, topo->rb_dihedrals[i], positions, forces,
                  (*energies)[ScalarStructure::DIHEDRAL], energies);
    }
  }

  template<class TBoundaryConditions>
  inline void RBDihedralSystemForce<TBoundaryConditions>::
  parallelEvaluate(const GenericTopology *topo, const Vector3DBlock *positions,
                   Vector3DBlock *forces, ScalarStructure *energies) {
    const TBoundaryConditions &boundary =
      (dynamic_cast<const SemiGenericTopology<TBoundaryConditions> &>(*topo)).
        boundaryConditions;

    unsigned int n = topo->rb_dihedrals.size();
    unsigned int count = numberOfBlocks(topo, positions);

    for (unsigned int i = 0; i < count; i++)
      if (Parallel::next()) {
        int to = (n * (i + 1)) / count;
        if (to > static_cast<int>(n))
          to = n;
        int from = (n * i) / count;
        for (int j = from; j < to; j++)
          calcRBTorsion(boundary, topo->rb_dihedrals[j], positions, forces,
                      (*energies)[ScalarStructure::DIHEDRAL], energies);
      }
  }

  template<class TBoundaryConditions>
  inline unsigned int RBDihedralSystemForce<TBoundaryConditions>::
  numberOfBlocks(const GenericTopology *topo, const Vector3DBlock *) {
    return std::min(Parallel::getAvailableNum(),
                    static_cast<int>(topo->rb_dihedrals.size()));
  }
}
#endif /* RBDihedralSystemForce_H */
