/* -*- c++ -*- */
#ifndef MTORSIONSYSTEMFORCE_H
#define MTORSIONSYSTEMFORCE_H

#include <protomol/force/system/SystemForce.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/topology/Torsion.h>
#include <protomol/type/Vector3DBlock.h>

namespace ProtoMol {
  //____ MTorsionSystemForce
  template<class TBoundaryConditions>
  class MTorsionSystemForce : public SystemForce {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class MTorsionSystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    void calcTorsion(const TBoundaryConditions &boundary,
                     const Torsion &currentTorsion,
                     const Vector3DBlock *positions,
                     Vector3DBlock *forces, Real &energy,
                     ScalarStructure *energies);
    Real calcTorsionEnergy(const TBoundaryConditions &boundary,
                           const Torsion &currentTorsion,
                           const Vector3DBlock *positions);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
  };

  //____ INLINES
  template<class TBoundaryConditions>
  inline void MTorsionSystemForce<TBoundaryConditions>::
  calcTorsion(const TBoundaryConditions &boundary, const Torsion &currTorsion,
              const Vector3DBlock *positions, Vector3DBlock *forces,
              Real &energy, ScalarStructure *energies) {

    int a1 = currTorsion.atom1;
    int a2 = currTorsion.atom2;
    int a3 = currTorsion.atom3;
    int a4 = currTorsion.atom4;

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
    for (int i = 0; i < currTorsion.multiplicity; i++)

      if (currTorsion.periodicity[i] > 0) {
        dpotdphi -= currTorsion.periodicity[i]
                    * currTorsion.forceConstant[i]
                    * sin(currTorsion.periodicity[i] * phi
                          + currTorsion.phaseShift[i]);

        // Add energy
        energy += currTorsion.forceConstant[i] *
          (1.0 + cos(currTorsion.periodicity[i] * phi +
                     currTorsion.phaseShift[i]));
      } else {
        Real diff = phi - currTorsion.phaseShift[i];

        if (diff < -M_PI)
          diff += 2 * M_PI;
        else if (diff > M_PI)
          diff -= 2 * M_PI;

        dpotdphi += 2.0 * currTorsion.forceConstant[i] * diff;

        // Add energy
        energy += currTorsion.forceConstant[i] * diff * diff;
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

  template<class TBoundaryConditions>
  inline Real MTorsionSystemForce<TBoundaryConditions>::
  calcTorsionEnergy(const TBoundaryConditions &boundary,
                    const Torsion &currTorsion,
                    const Vector3DBlock *positions) {
    int a1 = currTorsion.atom1;
    int a2 = currTorsion.atom2;
    int a3 = currTorsion.atom3;
    int a4 = currTorsion.atom4;

    Vector3D r12 = boundary.minimalDifference((*positions)[a2],
                                              (*positions)[a1]);
    Vector3D r23 = boundary.minimalDifference((*positions)[a3],
                                              (*positions)[a2]);
    Vector3D r34 = boundary.minimalDifference((*positions)[a4],
                                              (*positions)[a3]);

    // Cross product of r12 and r23, represents the plane shared by these two 
    // vectors
    Vector3D a = r12.cross(r23);
    // Cross product of r12 and r23, represents the plane shared by these two
    // vectors
    Vector3D b = r23.cross(r34);

    Vector3D c = r23.cross(a);

    // Calculate phi.
    Real cosPhi = a.dot(b) / (a.norm() * b.norm());
    Real sinPhi = c.dot(b) / (c.norm() * b.norm());
    Real phi = -atan2(sinPhi, cosPhi);

    // Calculate energy.

    Real energy = 0.0;

    for (int i = 0; i < currTorsion.multiplicity; i++)

      if (currTorsion.periodicity[i] > 0)

        // Add energy
        energy += currTorsion.forceConstant[i] *
          (1.0 + cos(currTorsion.periodicity[i] * phi +
                     currTorsion.phaseShift[i]));


      else {
        Real diff = phi - currTorsion.phaseShift[i];

        if (diff < -M_PI) diff += 2 * M_PI;
        else if (diff > M_PI) diff -= 2 * M_PI;

        // Add energy
        energy += currTorsion.forceConstant[i] * diff * diff;
      }

    return energy;
  }
}
#endif /* MTORSIONSYSTEMFORCE_H */
