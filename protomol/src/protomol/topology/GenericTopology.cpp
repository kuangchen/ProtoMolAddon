#include <protomol/topology/GenericTopology.h>

#include <protomol/base/Exception.h>

using namespace std;
using namespace ProtoMol;

//____ GenericTopology

const string GenericTopology::scope("Topology");
const string GenericTopology::keyword("Topology");

GenericTopology::GenericTopology() :
  exclude(ExclusionType::ONE4MODIFIED), coulombScalingFactor(1.0), time(0.0),
  min(Vector3D(Constant::MAXREAL, Constant::MAXREAL, Constant::MAXREAL)),
  max(Vector3D(-Constant::MINREAL, -Constant::MINREAL, -Constant ::MINREAL)),
  implicitSolvent(NONE), doSCPISM(0), forceFieldFlag(CHARMM), 
  minimalMolecularDistances(false), doGBSAOpenMM(0), obcType(0),
  dielecOffset(0), alphaObc(0), betaObc(0), gammaObc(0) {}

GenericTopology::GenericTopology(Real c, const ExclusionType &e) :
  exclude(e), coulombScalingFactor(c), time(0.0),
  min(Vector3D(Constant::MAXREAL, Constant::MAXREAL, Constant::MAXREAL)),
  max(Vector3D(-Constant::MINREAL, -Constant::MINREAL, -Constant::MINREAL)),
  implicitSolvent(NONE), doSCPISM(0), forceFieldFlag(CHARMM), 
  minimalMolecularDistances(false), doGBSAOpenMM(0), obcType(0),
  dielecOffset(0), alphaObc(0), betaObc(0), gammaObc(0) {}

GenericTopology *GenericTopology::make(const vector<Value> &values) const {
  assertParameters(values);

  return adjustAlias(doMake(values));
}

