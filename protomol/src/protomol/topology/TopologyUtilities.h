/* -*- c++ -*- */
#ifndef TOPOLOGYUTILITIES_H
#define TOPOLOGYUTILITIES_H

#include <protomol/type/Vector3DBlock.h>
#include <protomol/type/Vector3D.h>
#include <protomol/type/Matrix3By3.h>
#include <protomol/type/Real.h>
#include <protomol/topology/Torsion.h>
#include <protomol/type/Stack.h>
#include <protomol/topology/AngleInfo.h>
#include <protomol/topology/Bond.h>

#include <string>
#include <vector>
#include <set>



namespace ProtoMol {
  enum {X_AXIS, Y_AXIS, Z_AXIS};

  class GenericTopology;
  class ScalarStructure;

  //________________________________________randomVelocity
  /// Gaussian random distributed velocities
  void randomVelocity(Real temperature,
                      const GenericTopology *topology,
                      Vector3DBlock *velocities,
                      unsigned int seed = 1234);

  //________________________________________randomVelocity
  /// Gaussian random distributed velocities rescaled to the given interval
  /// with optional re-movements of linear and/or angular momentum
  void randomVelocity(Real temperatureFrom,
                      Real temperatureTo,
                      const GenericTopology *topology,
                      const Vector3DBlock *positions,
                      Vector3DBlock *velocities,
                      bool removeLinear = false,
                      bool removeAngular = false,
                      unsigned int seed = 1234);

  //________________________________________getAtomsBondedtoDihedral
  /// this function gets all the atoms bonded to ONE side of the dihedral

  void getAtomsBondedtoDihedral(const GenericTopology *topology,
                                std::set<int> *atomSet,
                                const int atomID,
                                const int inAtomID,
                                const int outAtomID,
                                const int exclAtomID);

  //________________________________________rotateDihedral
  /// this function rotates all the atoms bonded to ONE side of the dihedral

  void rotateDihedral(const GenericTopology *topology,
                      Vector3DBlock *positions,
                      const int dihedralID,
                      Real angle);

  void rotateDihedral(const GenericTopology *topology,
                      Vector3DBlock *positions,
                      Vector3DBlock *velocities,
                      const int dihedralID,
                      Real angle);

  //________________________________________kineticEnergy
  Real kineticEnergy(const GenericTopology *topology,
                     const Vector3DBlock *velocities);

  //________________________________________molecularKineticEnergy
  Real molecularKineticEnergy(const GenericTopology *topology,
                              const Vector3DBlock *velocitiies);

  //________________________________________kineticEnergyForAtomType
  enum waterOption {IGNORE_WATER, ONLY_WATER, ALL};
  Real kineticEnergyForAtomType(const GenericTopology *topology,
                                const Vector3DBlock *velocities,
                                int atomType,
                                waterOption option);
  Real kineticEnergyForAtomType(const GenericTopology *topology,
                                const Vector3DBlock *velocities,
                                int atomType,
                                waterOption option,
                                int &atomCount);

  //________________________________________kineticEnergyForWater
  Real kineticEnergyForWater(const GenericTopology *topology,
                             const Vector3DBlock *velocities);
  Real kineticEnergyForWater(const GenericTopology *topology,
                             const Vector3DBlock *velocities,
                             int &waterCount);

  //________________________________________kineticEnergyForNonWater
  Real kineticEnergyForNonWater(const GenericTopology *topology,
                                const Vector3DBlock *velocities);
  Real kineticEnergyForNonWater(const GenericTopology *topology,
                                const Vector3DBlock *velocities,
                                int &nonWaterCount);

  //________________________________________temperature
  Real temperature(const GenericTopology *topology,
                   const Vector3DBlock *velocities);

  //________________________________________temperature
  Real temperature(Real kineticEnergy, unsigned int degreesOfFreedom);

  //________________________________________temperatureForAtomType
  Real temperatureForAtomType(const GenericTopology *topology,
                              const Vector3DBlock *velocities,
                              int atomType,
                              waterOption option);

  //________________________________________temperatureForWater
  Real temperatureForWater(const GenericTopology *topology,
                           const Vector3DBlock *velocities);

  //________________________________________temperatureForNonWater
  Real temperatureForNonWater(const GenericTopology *topology,
                              const Vector3DBlock *velocities);

  //________________________________________getNonWaterAtoms
  int getNonWaterAtoms(const GenericTopology *topology);

  //________________________________________atomTypeToSymbolName
  std::string atomTypeToSymbolName(const std::string &type);

  //________________________________________linearMomentum
  Vector3D linearMomentum(const Vector3DBlock *velocities,
                          const GenericTopology *topo);

  //________________________________________linearMomentumSolute
  Vector3D linearMomentumSolute(const Vector3DBlock *velocities,
                                const GenericTopology *topo);

  //________________________________________centerOfMass
  Vector3D centerOfMass(const Vector3DBlock *positions,
                        const GenericTopology *topo);

  //________________________________________angularMomentum
  Vector3D angularMomentum(const Vector3DBlock *positions,
                           const Vector3DBlock *velocities,
                           const GenericTopology *topo);

  //________________________________________angularMomentum
  Vector3D angularMomentum(const Vector3DBlock *positions,
                           const Vector3DBlock *velocities,
                           const GenericTopology *topo,
                           const Vector3D &centerOfMass);

  //________________________________________angularMomentumSolute
  Vector3D angularMomentumSolute(const Vector3DBlock *positions,
                                 const Vector3DBlock *velocities,
                                 const GenericTopology *topo,
                                 const Vector3D &centerOfMass);

  //________________________________________inertiaMomentum
  Matrix3By3 inertiaMomentum(const Vector3DBlock *positions,
                             const GenericTopology *topo,
                             const Vector3D &centerOfMass);

  //________________________________________inertiaMomentumSolute
  Matrix3By3 inertiaMomentumSolute(const Vector3DBlock *positions,
                                   const GenericTopology *topo,
                                   const Vector3D &centerOfMass);

  //________________________________________removeLinearMomentum
  Vector3D removeLinearMomentum(Vector3DBlock *velocities,
                                const GenericTopology *topo);

  //________________________________________removeAngularMomentum
  Vector3D removeAngularMomentum(const Vector3DBlock *positions,
                                 Vector3DBlock *velocities,
                                 const GenericTopology *topo);

  //________________________________________velocityVirial
  ScalarStructure velocityVirial(const GenericTopology *topology,
                                 const Vector3DBlock *velocities);

  //________________________________________addVelocityVirial
  void addVelocityVirial(ScalarStructure *energies,
                         const GenericTopology *topology,
                         const Vector3DBlock *velocities);

  //________________________________________computePressure
  Real computePressure(const GenericTopology *topology,
                       const Vector3DBlock *positions,
                       const Vector3DBlock *velocities,
                       const ScalarStructure *energies);

  //________________________________________computePressure
  Real computePressure(const ScalarStructure *energies,
                       Real volume,
                       Real kineticEnergy);

  //________________________________________computeMolecularPressure
  Real computeMolecularPressure(const ScalarStructure *energies,
                                Real volume,
                                Real kineticEnergy);

  //________________________________________molecularMomentum
  Vector3D molecularMomentum(const std::vector<int> &,
                             const Vector3DBlock *,
                             const GenericTopology *);

  //________________________________________molecularCenterOfMass
  Vector3D molecularCenterOfMass(const std::vector<int> &,
                                 const Vector3DBlock *,
                                 const GenericTopology *);

  //________________________________________buildMolecularCenterOfMass
  void buildMolecularCenterOfMass(const Vector3DBlock *positions,
                                  GenericTopology *topo);

  //________________________________________buildMolecularMomentum
  void buildMolecularMomentum(const Vector3DBlock *velocities,
                              GenericTopology *topo);


  //________________________________________buildRattleShakeBondConstraintList
  void buildRattleShakeBondConstraintList(
    GenericTopology *topo,
    std::vector<Bond::Constraint> &
    bondConstraints, bool all);

  void build_angle_list(const GenericTopology *topo,
                        const unsigned int atomID,
                        const unsigned int inAtomID,
                        const unsigned int outAtomID,
                        const unsigned int exclAtomID,
                        Real rotAngle,
                        std::vector<AngleInfo> *angles);
  void set_angles(Stack<unsigned int> *nodeStack,
                  std::vector<AngleInfo> *angles,
                  bool lastIsInnerAtom,
                  Real wholeAngle);
  void general_rotation(unsigned int innerAtom1,
                        unsigned int innerAtom2,
                        Vector3DBlock *positions,
                        std::vector<AngleInfo> *angles);
  void general_rotation(unsigned int innerAtom1,
                        unsigned int innerAtom2,
                        Vector3DBlock *positions,
                        Vector3DBlock *velocities,
                        std::vector<AngleInfo> *angles);
  Real computePhiDihedral(const GenericTopology *topo,
                          const Vector3DBlock *positions,
                          int index);
  Real computePhiDihedralEnergy(const GenericTopology *topo,
                                int index,
                                Real phi);
}

#endif // TOPOLOGYUTILITIES_H

