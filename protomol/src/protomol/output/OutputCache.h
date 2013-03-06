/*  -*- c++ -*-  */
#ifndef PROTOMOL_OUTPUT_CACHE_H
#define PROTOMOL_OUTPUT_CACHE_H

#include <protomol/type/Real.h>
#include <protomol/type/PDB.h>
#include <protomol/type/PAR.h>
#include <protomol/type/PSF.h>
#include <protomol/type/Vector3D.h>

namespace ProtoMol {
  class ProtoMolApp;
  class Output;
  class Configuration;
  class GenericTopology;
  class ScalarStructure;
  class Vector3DBlock;
  class OutputFactory;
  class Integrator;
  struct PDB;
  struct Atom;

  /**
     OutputCache caches all kind of values, which may be needed
     by Output objects and simplifies the access to values of interest.
     Add new cached values, if needed ..
     There are some (feature) values, which will only change when th
     Topology chang
   */
  class OutputCache  {
    const ProtoMolApp *app;

    Vector3DBlock *initialPositions;
    mutable Vector3DBlock * minimalPositions;

    //  Additional dat
    std::vector<PDB::Atom> atoms;
    std::vector<Real> remRates;
    const Real *replicaHistory;
    PSF psf;
    PAR par;

    mutable bool cachedKE;
    mutable Real kE;
    mutable Real t;

    mutable bool cachedPE;
    mutable Real pE;

    mutable bool cachedV;
    mutable Real v;

    mutable bool cachedP;
    mutable Real p;
    mutable bool cachedMolP;
    mutable Real molP;

    mutable bool cachedLinearMomentum;
    mutable Vector3D linearMomentum;

    mutable bool cachedAngularMomentum;
    mutable Vector3D angularMomentum;

    mutable bool cachedCenterOfMass;
    mutable Vector3D centerOfMass;

    mutable bool cachedDiffusion;
    mutable Real diffusion;

    mutable bool cachedDensity;
    mutable Real density;

    mutable bool cachedMass;
    mutable Real mass;

    mutable int cachedDihedralPhi;
    mutable Real dihedralPhi;

    mutable bool cachedDihedralPhis;
    mutable std::vector<Real> *dihedralPhis;

    mutable bool cachedBrentMaxima;
    mutable std::vector<std::vector<Real> > *brentMaxima;

    mutable bool cachedMolT;
    mutable Real molT;

    mutable bool cachedMolKE;
    mutable Real molKE;

    mutable bool cachedMinimalPositions;

    bool restore;

  public:
    OutputCache();
    ~OutputCache();

    void initialize(const ProtoMolApp *app);

    //  Methods to add additional data for output object
    void add(const std::vector<PDB::Atom> &pdbAtoms) {atoms = pdbAtoms;}
    void add(const PSF &psf) {this->psf = psf;}
    void add(const PAR &par) {this->par = par;}
    void add(const std::vector<Real> &remRates) {this->remRates = remRates;}
    void add(const Real *replicaHistory)
    {this->replicaHistory = replicaHistory;}

    Real getTotalEnergy() const;
    Real getPotentialEnergy() const;
    Real getKineticEnergy() const;
    Real getTemperature() const;
    Real getVolume() const;
    Real getTime() const;
    Real getPressure() const;
    Real getMolecularPressure() const;
    Real getMolecularTemperature() const;
    Real getMolecularKineticEnergy() const;
    Vector3D getLinearMomentum() const;
    Vector3D getAngularMomentum() const;
    Vector3D getCenterOfMass() const;
    Real getDiffusion() const;
    Real getDensity() const;
    Real getMass() const;
    Real getDihedralPhi(int index) const;
    Real getBrent(Real ax, Real bx, Real cx, Real tol, Real &xmin, int dihindex,
                  bool max) const;
    std::vector<Real> getDihedralPhis(std::vector<int> ) const;
    std::vector<std::vector<Real> > getBrentMaxima(std::vector<int> ,
                                                   bool) const;
    const Vector3DBlock *getMinimalPositions() const;
    const std::vector<Real> &getREMRates() const {return remRates;}
    const Real *getReplicaHistory() const {return replicaHistory;}

    const std::vector<PDB::Atom> &getAtoms() const {return atoms;}
    const PSF &getPSF() const {return psf;}
    const PAR &getPAR() const {return par;}

    // / To be called before every run() or finialize()
    void uncache() const;

    void setRestore() {restore = true;}
    void clearRestore() {restore = false;}
    bool getRestore() const {return restore;}

    const ProtoMolApp *getApp() const {return app;}
  };
}
#endif //  PROTOMOL_OUTPUT_CACHE_H
