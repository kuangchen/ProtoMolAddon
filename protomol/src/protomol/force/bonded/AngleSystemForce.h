/* -*- c++ -*- */
#ifndef ANGLESYSTEMFORCE_H
#define ANGLESYSTEMFORCE_H

#include <protomol/force/system/SystemForce.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/parallel/Parallel.h>

namespace ProtoMol {
  //____ AngleSystemForce
  template<class TBoundaryConditions>
  class AngleSystemForce : public SystemForce {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual ~AngleSystemForce() {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class AngleSystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    void calcAngle(const TBoundaryConditions &boundary,
                   const Angle &currentAngle, const Vector3DBlock *positions,
                   Vector3DBlock *forces,
                   ScalarStructure *energies);

    Real calcAngleEnergy(const TBoundaryConditions &boundary,
                         const Angle &currentAngle,
                         const Vector3DBlock *positions);

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
    virtual std::string getKeyword() const {return "Angle";}

    virtual unsigned int numberOfBlocks(const GenericTopology *topo,
                                        const Vector3DBlock *pos);

  private:
    virtual Force *doMake(const std::vector<Value> &) const {
      return new AngleSystemForce();
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const {return getKeyword();}
    virtual void getParameters(std::vector<Parameter> &) const {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
  };

  //____ INLINES

  template<class TBoundaryConditions>
  inline void AngleSystemForce<TBoundaryConditions>::evaluate(
    const GenericTopology *topo, const Vector3DBlock *positions,
    Vector3DBlock *forces, ScalarStructure *energies) {

    const TBoundaryConditions &boundary =
      ((SemiGenericTopology<TBoundaryConditions> &)(*topo)).
        boundaryConditions;
    for (unsigned int i = 0; i < topo->angles.size(); i++)
      calcAngle(boundary, topo->angles[i], positions, forces, energies);
  }

  template<class TBoundaryConditions>
  inline void AngleSystemForce<TBoundaryConditions>::calcAngle(
    const TBoundaryConditions &boundary, const Angle &currentAngle,
    const Vector3DBlock *positions, Vector3DBlock *forces,
    ScalarStructure *energies) {
    int a1 = currentAngle.atom1;
    int a2 = currentAngle.atom2;
    int a3 = currentAngle.atom3;

    Real restAngle = currentAngle.restAngle;
    Real forceConstant = currentAngle.forceConstant;
    Real ureyBradleyConstant = currentAngle.ureyBradleyConstant;
    Real ureyBradleyRestLength = currentAngle.ureyBradleyRestLength;

    Vector3D atom1((*positions)[a1]);
    Vector3D atom2((*positions)[a2]);
    Vector3D atom3((*positions)[a3]);

    Vector3D r12(boundary.minimalDifference(atom2, atom1));
    Vector3D r32(boundary.minimalDifference(atom2, atom3));
    Vector3D r13(boundary.minimalDifference(atom3, atom1));
    Real d12 = r12.norm();                     // Distance between atom 1 and 2.
    Real d32 = r32.norm();                     // Distance between atom 3 and 2.
    Real d13 = r13.norm();                     // Distance between atom 1 and 3.

    // Calculate theta.
    Real theta = atan2((r12.cross(r32)).norm(), r12.dot(r32));
    Real sinTheta = sin(theta);
    Real cosTheta = cos(theta);

    // Calculate dpot/dtheta
    Real dpotdtheta = 2.0 * forceConstant * (theta - restAngle);

    // Calculate dr/dx, dr/dy, dr/dz.
    Vector3D dr12(r12 / d12);
    Vector3D dr32(r32 / d32);
    Vector3D dr13(r13 / d13);
    // Calulate dtheta/dx, dtheta/dy, dtheta/dz.
    Vector3D dtheta1((dr12 * cosTheta - dr32) / (sinTheta * d12));     // atom1
    Vector3D dtheta3((dr32 * cosTheta - dr12) / (sinTheta * d32));     // atom3

    // Calculate Urey Bradley force.
    Vector3D ureyBradleyforce1(
      dr13 * (2.0 * ureyBradleyConstant * (d13 - ureyBradleyRestLength)));
    Vector3D ureyBradleyforce3(-ureyBradleyforce1);

    // Calculate force on atom1 due to atom 2 and 3.
    Vector3D force1(-dtheta1 * dpotdtheta - ureyBradleyforce1);

    // Calculate force on atom3 due to atom 1 and 2.
    Vector3D force3(-dtheta3 * dpotdtheta - ureyBradleyforce3);

    // Calculate force on atom2 due to atom 1 and 3.
    Vector3D force2(-force1 - force3);

    // Add to the total force.
    (*forces)[a1] += force1;
    (*forces)[a2] += force2;
    (*forces)[a3] += force3;

    // Calculate Energy.
    Real eHarmonic = forceConstant * (theta - restAngle) * (theta - restAngle);
    Real eUreyBradley = ureyBradleyConstant *
                        (d13 -
                         ureyBradleyRestLength) * (d13 - ureyBradleyRestLength);
    // Add Energy
    (*energies)[ScalarStructure::ANGLE] += eHarmonic + eUreyBradley;

    // Add virial
    if (energies->virial()) {
      Real xy = force1.c[0] * r12.c[1] + force3.c[0] * r32.c[1];
      Real xz = force1.c[0] * r12.c[2] + force3.c[0] * r32.c[2];
      Real yz = force1.c[1] * r12.c[2] + force3.c[1] * r32.c[2];
      (*energies)[ScalarStructure::VIRIALXX] += force1.c[0] * r12.c[0] + force3.c[0] *
                                                r32.c[0];
      (*energies)[ScalarStructure::VIRIALXY] += xy;
      (*energies)[ScalarStructure::VIRIALXZ] += xz;
      (*energies)[ScalarStructure::VIRIALYX] += xy;
      (*energies)[ScalarStructure::VIRIALYY] += force1.c[1] * r12.c[1] + force3.c[1] *
                                                r32.c[1];
      (*energies)[ScalarStructure::VIRIALYZ] += yz;
      (*energies)[ScalarStructure::VIRIALZX] += xz;
      (*energies)[ScalarStructure::VIRIALZY] += yz;
      (*energies)[ScalarStructure::VIRIALZZ] += force1.c[2] * r12.c[2] + force3.c[2] *
                                                r32.c[2];
    }
  }

  template<class TBoundaryConditions>
  inline Real AngleSystemForce<TBoundaryConditions>::calcAngleEnergy(
    const TBoundaryConditions &boundary, const Angle &currentAngle,
    const Vector3DBlock *positions) {
    int a1, a2, a3;
    Real restAngle, forceConstant, ureyBradleyConstant, ureyBradleyRestLength;
    Vector3D atom1, atom2, atom3, r12, r32, r13;
    Real d13, theta, eHarmonic, eUreyBradley;

    a1 = currentAngle.atom1;
    a2 = currentAngle.atom2;
    a3 = currentAngle.atom3;
    restAngle = currentAngle.restAngle;
    forceConstant = currentAngle.forceConstant;
    ureyBradleyConstant = currentAngle.ureyBradleyConstant;
    ureyBradleyRestLength = currentAngle.ureyBradleyRestLength;

    atom1 = (*positions)[a1];
    atom2 = (*positions)[a2];
    atom3 = (*positions)[a3];

    r12 = boundary.minimalDifference(atom2, atom1);
    r32 = boundary.minimalDifference(atom2, atom3);
    r13 = boundary.minimalDifference(atom3, atom1);
    d13 = r13.norm();                        // Distance between atom 1 and 3.

    // Calculate theta.
    theta = atan2((r12.cross(r32)).norm(), r12.dot(r32));

    // Calculate Energy.
    eHarmonic = forceConstant * (theta - restAngle) * (theta - restAngle);
    eUreyBradley = ureyBradleyConstant * (d13 - ureyBradleyRestLength)
                   * (d13 - ureyBradleyRestLength);

    return eHarmonic + eUreyBradley;
  }

  template<class TBoundaryConditions>
  inline void AngleSystemForce<TBoundaryConditions>::parallelEvaluate(
    const GenericTopology *topo, const Vector3DBlock *positions,
    Vector3DBlock *forces, ScalarStructure *energies) {
    const TBoundaryConditions &boundary =
      (dynamic_cast<const SemiGenericTopology<TBoundaryConditions> &>(*topo)).
        boundaryConditions;

    unsigned int n = topo->angles.size();
    unsigned int count = numberOfBlocks(topo, positions);

    for (unsigned int i = 0; i < count; i++)
      if (Parallel::next()) {
        int to = (n * (i + 1)) / count;
        if (to > static_cast<int>(n))
          to = n;
        int from = (n * i) / count;
        for (int j = from; j < to; j++)
          calcAngle(boundary, topo->angles[j], positions, forces, energies);
      }
  }

  template<class TBoundaryConditions>
  inline unsigned int AngleSystemForce<TBoundaryConditions>::numberOfBlocks(
    const GenericTopology *topo, const Vector3DBlock *) {
    return std::min(Parallel::getAvailableNum(),
                    static_cast<int>(topo->angles.size()));
  }
}
#endif /* ANGLESYSTEMFORCE_H */
