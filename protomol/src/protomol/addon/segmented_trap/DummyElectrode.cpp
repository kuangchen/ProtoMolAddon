#include <protomol/addon/segmented_trap/DummyElectrode.h>

using namespace ProtoMolAddon::SegmentedTrap;

DummyElectrode::DummyElectrode() {}

double DummyElectrode::GetPotential(const Vector3D &pos, double t) const {
  return 0;
}
