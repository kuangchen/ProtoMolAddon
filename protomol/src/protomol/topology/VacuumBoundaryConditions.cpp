#include <protomol/topology/VacuumBoundaryConditions.h>

using namespace std;
using namespace ProtoMol;
//____ VacuumBoundaryConditions

const string VacuumBoundaryConditions::keyword("Vacuum");

Vector3D VacuumBoundaryConditions::getMin() const {
  return Vector3D(-1.0 * Constant::REAL_INFINITY,
    -1.0 * Constant::REAL_INFINITY,
    -1.0 * Constant::REAL_INFINITY);
}

Vector3D VacuumBoundaryConditions::getMax() const {
  return Vector3D(Constant::REAL_INFINITY,
    Constant::REAL_INFINITY,
    Constant::REAL_INFINITY);
}

Real VacuumBoundaryConditions::getVolume() const {
  return Constant::REAL_INFINITY;
}

Vector3D VacuumBoundaryConditions::e1()     const {
  return Vector3D(0, 0, 0);
}

Vector3D VacuumBoundaryConditions::e2()     const {
  return Vector3D(0, 0, 0);
}

Vector3D VacuumBoundaryConditions::e3()     const {
  return Vector3D(0, 0, 0);
}

Vector3D VacuumBoundaryConditions::e1r()    const {
  return Vector3D(0, 0, 0);
}

Vector3D VacuumBoundaryConditions::e2r()    const {
  return Vector3D(0, 0, 0);
}

Vector3D VacuumBoundaryConditions::e3r()    const {
  return Vector3D(0, 0, 0);
}

Vector3D VacuumBoundaryConditions::origin() const {
  return Vector3D(0, 0, 0);
}

vector<Vector3D> VacuumBoundaryConditions::buildLatticeVectors(Real) const {
  return vector<Vector3D>();
}

VacuumBoundaryConditions VacuumBoundaryConditions::make(vector<Value> ) {
  return VacuumBoundaryConditions();
}

