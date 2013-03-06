#include <protomol/force/CompareForce.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/base/MathUtilities.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/parallel/Parallel.h>

using namespace std;
using namespace ProtoMol::Report;

using namespace ProtoMol;
//____ CompareForce

const string CompareForce::keyword("compare");
unsigned int CompareForce::myCounter = 0;

CompareForce::CompareForce(Force *actualForce,
                           CompareForce *compareForce) :
  myForces(new Vector3DBlock), myEnergies(new ScalarStructure) {
  if (actualForce == NULL)
    report << error << "Actual force is a zero pointer." << endr;
  myActualForce = actualForce;
  myCompareForce = compareForce;
  if (compareForce != NULL) {
    myCompareForcename = compareForce->getId() + "." + toString(
      myCompareForce->getIdNumber());
    myForcename = actualForce->getId() + "." + toString(myCounter);
  }
  myIdNumber = myCounter;
  myCounter++;
}

CompareForce::~CompareForce() {
  unsigned int n = myErrors.size();
  if (n > 0) {
    CompareError avg;
    for (unsigned int i = 0; i < n; ++i) {
      avg.absF2 += myErrors[i].absF2 / n;
      avg.rFavg += myErrors[i].rFavg / n;
      avg.rFmax += myErrors[i].rFmax / n;
      avg.rPE += myErrors[i].rPE / n;
    }

    CompareError stddev;
    if (n > 1) {
      for (unsigned int i = 0; i < n; ++i) {
        stddev.absF2 += power<2>(myErrors[i].absF2 - avg.absF2);
        stddev.rFavg += power<2>(myErrors[i].rFavg - avg.rFavg);
        stddev.rFmax += power<2>(myErrors[i].rFmax - avg.rFmax);
        stddev.rPE += power<2>(myErrors[i].rPE - avg.rPE);
      }

      Real m = n - 1.0;
      stddev.absF2 = sqrt(stddev.absF2 / m);
      stddev.rFavg = sqrt(stddev.rFavg / m);
      stddev.rFmax = sqrt(stddev.rFmax / m);
      stddev.rPE = sqrt(stddev.rPE / m);
    }
    report << plain << "Comparing " << toString(myIdNumber / 2) << " :"
           << " absF2 " << toString(avg.absF2) << " (" << toString(
      stddev.absF2) << ")"
           << ", rFavg " << toString(avg.rFavg) << " (" << toString(
      stddev.rFavg) << ")"
           << ", rFmax " << toString(avg.rFmax) << " (" << toString(
      stddev.rFmax) << ")"
           << ", rPE " << toString(avg.rPE) << " (" <<
    toString(stddev.rPE) << ")"
           << ", n=" << n
           << endr;
  }

  if (myActualForce != NULL) delete myActualForce;
  if (myForces != NULL) delete myForces;
  if (myEnergies != NULL) delete myEnergies;
}

void CompareForce::preprocess(unsigned int numAtoms) {
  myForces->zero(numAtoms);
  myEnergies->clear();
}

void CompareForce::postprocess(const GenericTopology *topo,
                               Vector3DBlock *forces,
                               ScalarStructure *energies) {
  if (myCompareForce == NULL) {
    // Add the contributions to the right place
    forces->intoAdd(*myForces);
    energies->intoAdd(*myEnergies);

  } else if (!Parallel::isParallel()) {
       // Compute the error of your choice.
    const ScalarStructure *comparedEnergies = myCompareForce->getEnergies();
    const Vector3DBlock *comparedForces = myCompareForce->getForces();
    Real difftotal = 0.0;
    Real total = 0.0;
    Real errmax = 0.0;
    Real abserrmax = 0.0;
    for (unsigned int i = 0; i < forces->size(); i++) {
      Real rm = 1.0 / topo->atoms[i].scaledMass;
      Real diff = ((*myForces)[i] - (*comparedForces)[i]).normSquared();
      Real errfmag = sqrt(diff * rm);
      Real fmag = sqrt((*myForces)[i].normSquared() * rm);
      if (errfmag > errmax)
        errmax = errfmag;
      difftotal += errfmag;
      total += fmag;
      if (diff > abserrmax)
        abserrmax = diff;
    }

    CompareError theError(sqrt(abserrmax),
                          (total != 0.0 ? difftotal / total : 0.0),
                          (total != 0.0 ? forces->size() * errmax /
                           total : 0.0),
                          (myEnergies->potentialEnergy() != 0.0 ?
                           fabs(
                             (comparedEnergies->potentialEnergy() -
                              myEnergies->potentialEnergy()) /
                             myEnergies->potentialEnergy()) :
                           0.0));

    myErrors.push_back(theError);

    report << plain
           << "Comparing " << toString(myIdNumber / 2) << " :"
           << " absF2 " << toString(theError.absF2)
           << ", rFavg " << toString(theError.rFavg)
           << ", rFmax " << toString(theError.rFmax)
           << ", rPE " << toString(theError.rPE)
           << endr;
  }
}

void CompareForce::getParameters(vector<Parameter> &parameters) const {
  myActualForce->getParameters(parameters);
}

string CompareForce::getIdNoAlias() const {
  return string(CompareForce::keyword + " " + myActualForce->getIdNoAlias());
}

Force *CompareForce::doMake(const vector<Value> &values) const {
  return myActualForce->make(values);
}

unsigned int CompareForce::numberOfBlocks(const GenericTopology *topo,
                                          const Vector3DBlock *positions) {
  return myActualForce->numberOfBlocks(topo, positions);
}

void CompareForce::uncache() {
  myActualForce->uncache();
}

void CompareForce::doSetParameters(vector<Value> values) {
  myActualForce->setParameters(values);
}

