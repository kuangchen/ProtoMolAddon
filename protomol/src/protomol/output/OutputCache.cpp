#include <protomol/output/OutputCache.h>
#include <protomol/output/Output.h>
#include <protomol/config/Configuration.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/integrator/Integrator.h>

#include <protomol/base/Report.h>
#include <protomol/base/MathUtilities.h>
#include <protomol/base/SystemUtilities.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>

#include <protomol/ProtoMolApp.h>

#include <protomol/base/Zap.h>
#include <protomol/base/Exception.h>

#include <limits>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

OutputCache::OutputCache() :
  initialPositions(new Vector3DBlock()),
  minimalPositions(new Vector3DBlock()),
  dihedralPhi(Constant::REAL_NAN),
  dihedralPhis(new vector<Real>()),
  brentMaxima(new vector<vector<Real> >()) {
  uncache();
}


OutputCache::~OutputCache() {
  zap(initialPositions);
  zap(minimalPositions);
  zap(dihedralPhis);
  zap(brentMaxima);
}


void OutputCache::initialize(const ProtoMolApp *app) {
  this->app = app;
  *initialPositions = app->positions;
}


Real OutputCache::getTotalEnergy() const {
  return getPotentialEnergy() + getKineticEnergy();
}


Real OutputCache::getPotentialEnergy() const {
  if (!cachedPE) {
    pE = app->energies.potentialEnergy();
    cachedPE = true;
  }
  return pE;
}


Real OutputCache::getKineticEnergy() const {
  if (!cachedKE) {
    kE = ProtoMol::kineticEnergy(app->topology, &app->velocities);
    t = ProtoMol::temperature(kE, app->topology->degreesOfFreedom);
    cachedKE = true;
  }
  return kE;
}


Real OutputCache::getTemperature() const {
  if (!cachedKE) {
    kE = ProtoMol::kineticEnergy(app->topology, &app->velocities);
    t = ProtoMol::temperature(kE, app->topology->degreesOfFreedom);
    cachedKE = true;
  }
  return t;
}


Real OutputCache::getMolecularTemperature() const {
  if (!cachedMolT) {
    molT = ProtoMol::temperature(getMolecularKineticEnergy(), 3 *
                                 app->topology->molecules.size());
    cachedMolT = true;
  }

  return molT;
}


Real OutputCache::getMolecularKineticEnergy() const {
  if (!cachedMolKE) {
    molKE = ProtoMol::molecularKineticEnergy(app->topology, &app->velocities);
    cachedMolKE = true;
  }

  return molKE;
}


Real OutputCache::getPressure() const {
  if (!cachedP) {
    if (!app->energies.virial())
      p = 0.0;
    else if (getVolume() > 0.0)
      p = ProtoMol::computePressure(&app->energies, getVolume(),
                                    getKineticEnergy());
    else p = Constant::REAL_INFINITY;
    cachedP = true;
  }

  return p;
}


Real OutputCache::getMolecularPressure() const {
  if (!cachedMolP) {
    if (!app->energies.molecularVirial())
      molP = 0.0;
    else if (getVolume() > 0.0)
      molP = ProtoMol::computeMolecularPressure(&app->energies,
                                                getVolume(),
                                                getMolecularKineticEnergy());
    else molP = Constant::REAL_INFINITY;
    cachedMolP = true;
  }

  return molP;
}


Real OutputCache::getVolume() const {
  if (!cachedV) {
    v = app->topology->getVolume(app->positions);
    cachedV = true;
  }
  return v;
}


Vector3D OutputCache::getLinearMomentum() const {
  if (!cachedLinearMomentum) {
    linearMomentum =
      ProtoMol::linearMomentum(&app->velocities, app->topology);
    cachedLinearMomentum = true;
  }
  return linearMomentum;
}

Vector3D OutputCache::getAngularMomentum() const {
  if (!cachedAngularMomentum) {
    angularMomentum =
      ProtoMol::angularMomentum(&app->positions, &app->velocities,
                                app->topology, getCenterOfMass());
    cachedAngularMomentum = true;
  }

  return angularMomentum;
}


Vector3D OutputCache::getCenterOfMass() const {
  if (!cachedCenterOfMass) {
    centerOfMass = ProtoMol::centerOfMass(&app->positions, app->topology);
    cachedCenterOfMass = true;
  }

  return centerOfMass;
}


Real OutputCache::getDiffusion() const {
  if (!cachedDiffusion) {
    diffusion = 0.0;
    unsigned int numberOfAtoms = app->positions.size();
    for (unsigned int i = 0; i < numberOfAtoms; i++)
      diffusion +=
        (app->positions[i] - (*initialPositions)[i]).normSquared();

    diffusion /= (6.0 * numberOfAtoms);
    cachedDiffusion = true;
  }

  return diffusion;
}


Real OutputCache::getDensity() const {
  if (!cachedDensity) {
    density =
      (getVolume() > 0 ?
       (getMass() / getVolume() * Constant::SI::AMU *
        power<3>(Constant::SI::LENGTH_AA) *
        1e-3) : Constant::REAL_NAN);

    cachedDensity = true;
  }

  return density;
}


Real OutputCache::getMass() const {
  if (!cachedMass) {
    mass = 0.0;
    unsigned int numberOfAtoms = app->positions.size();
    for (unsigned int i = 0; i < numberOfAtoms; i++)
      mass += app->topology->atoms[i].scaledMass;

    cachedMass = true;
  }

  return mass;
}


Real OutputCache::getTime() const {
  return app->topology->time;
}


const Vector3DBlock *OutputCache::getMinimalPositions() const {
  if (!cachedMinimalPositions) {
    *minimalPositions = app->positions;
    (const_cast<GenericTopology *>(app->topology))->
      minimalImage(*minimalPositions);
  }
  cachedMinimalPositions = true;

  return minimalPositions;
}


Real OutputCache::getDihedralPhi(int index) const {
  if (index < 0 || index >= static_cast<int>(app->topology->dihedrals.size()))
    index = -1;

  if (index < 0) {
    dihedralPhi = Constant::REAL_NAN;
    cachedDihedralPhi = index;
    return dihedralPhi;
  }

  dihedralPhi = computePhiDihedral(app->topology, &app->positions, index);

  return dihedralPhi;
}


vector<Real> OutputCache::getDihedralPhis(vector<int> dihedralset) const {
  if (!cachedDihedralPhis) {
    dihedralPhis->resize(dihedralset.size());
    for (unsigned int i = 0; i < dihedralset.size(); ++i)
      (*dihedralPhis)[i] =
        computePhiDihedral(app->topology, &app->positions, dihedralset[i]);

    //  different functions require different dihedralset
    cachedDihedralPhis = false;
  }

  return *dihedralPhis;
}


//  Brent's Maxima function goes here to use topology for dihedral well
//  calculation
vector<vector<Real> > OutputCache::getBrentMaxima(vector<int> dihedralset,
                                                  bool max) const {
  if (!cachedBrentMaxima) {
    // The Brent algorithm gets the maxima if maxmin = -1
    int maxmin = -1;
    if (!max)
      maxmin = 1;

    brentMaxima->clear();
    brentMaxima->resize(dihedralset.size());

    for (unsigned int i = 0; i < dihedralset.size(); ++i) {
      // note the function evaluates one step past 2 pi
      for (unsigned int j = 0; j <= 99; j++) {
        Real lradangle = (M_PI * 2 / 100 * j);
        Real radangle = (M_PI * 2 / 100 * (j + 1));
        Real rradangle = (M_PI * 2 / 100 * (j + 2));

        Real valLangle = maxmin *
          computePhiDihedralEnergy(app->topology, dihedralset[i], lradangle);

        Real valRangle = maxmin *
          computePhiDihedralEnergy(app->topology, dihedralset[i], rradangle);

        Real valAngle = maxmin *
          computePhiDihedralEnergy(app->topology, dihedralset[i], radangle);

        Real xmax = 0.0;
        Real tol = 0.01;
        if ((valLangle > valAngle) && (valRangle > valAngle)) {
          getBrent(lradangle, radangle, rradangle, tol, xmax, dihedralset[i],
                   max);

          ((*brentMaxima)[i]).push_back(xmax);
        }
      }

      //  Throws Warning if no maxima were found
      if (((*brentMaxima)[i]).size() == 0) {
        report << warning
               << "No dihedral maxima found for dihedral index: "
               << dihedralset[i] << " Check dihedral energy equation"
               << endr;
        ((*brentMaxima)[i]).push_back(0.0);
      }
    }

    // temp hack that allows multiple calls but defeats the purpose of
    //  the app.outputCache... please fix me!
    // hedBrentMaxima= true;
  }

  return *brentMaxima;
}


// BRENT FUNCTION
Real OutputCache::getBrent(Real ax, Real bx, Real cx, Real tol, Real &xmin,
                           int dihindex, bool max) const {
  const int ITMAX = 100;
  const Real CGOLD = 0.3819660;
  const Real ZEPS = numeric_limits<Real>::epsilon() * 1.0e-3;
  Real a, b, d = 0.0, etemp, fu, fv, fw, fx;
  Real p, q, r, tol1, tol2, u, v, w, x, xm;
  Real e = 0.0;

  // The Brent algorithm gets the maxima if maxmin = -1
  int maxmin = -1;
  if (!max) maxmin = 1;

  a = (ax < cx ? ax : cx);
  b = (ax > cx ? ax : cx);
  x = w = v = bx;
  fw = fv = fx = maxmin * computePhiDihedralEnergy(app->topology, dihindex, x);
  for (int iter = 0; iter < ITMAX; iter++) {
    xm = 0.5 * (a + b);
    tol2 = 2.0 * (tol1 = tol * fabs(x) + ZEPS);
    if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) {
      xmin = x;
      return fx;
    }

    if (fabs(e) > tol1) {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q < 0.0) p = -p;
      q = fabs(q);
      etemp = e;
      e = d;

      if (fabs(p) >=
          fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
        d = CGOLD * (e = (x >= xm ? a - x : b - x));

      else {
        d = p / q;
        u = x + d;
        if (u - a < tol2 || b - u < tol2)
          d = sign(tol1, xm - x);
      }

    } else d = CGOLD * (e = (x >= xm ? a - x : b - x));
    u = (fabs(d) >= tol1 ? x + d : x + sign(tol1, d));
    fu = maxmin * computePhiDihedralEnergy(app->topology, dihindex, u);

    if (fu <= fx) {
      if (u >= x) a = x;else b = x;
      shift(v, w, x, u);
      shift(fv, fw, fx, fu);

    } else {
      if (u < x) a = u;else b = u;
      if (fu <= fw || w == x) {
        v = w;
        w = u;
        fv = fw;
        fw = fu;

      } else if (fu <= fv || v == x || v == w) {
        v = u;
        fv = fu;
      }
    }
  }

  // too many iterations in brent
  xmin = x;

  return fx;
}


void OutputCache::uncache() const {
  cachedKE = false;
  cachedPE = false;
  cachedV = false;
  cachedP = false;
  cachedMolP = false;
  cachedLinearMomentum = false;
  cachedAngularMomentum = false;
  cachedCenterOfMass = false;
  cachedDiffusion = false;
  cachedDensity = false;
  cachedMass = false;
  cachedDihedralPhis = false;
  cachedDihedralPhi = -1;
  cachedBrentMaxima = false;
  cachedMolT = false;
  cachedMolKE = false;
  cachedMinimalPositions = false;
}
