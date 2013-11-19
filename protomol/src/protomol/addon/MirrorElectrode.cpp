#include <protomol/addon/MirrorElectrode.h>

using namespace ProtoMol;
using namespace ProtoMolAddon;

MirrorElectrode::MirrorElectrode(const string& label, const RealElectrode* parent, const std::string& reflection):
  AbstractElectrode(label),
  parent(parent),
  reflection(reflection) {
}
           
double MirrorElectrode::GetNNPotential(const Vector3D& pos, const array<int, 3>& offset) const {
  //cout << "label = " << label << "pos = " << pos << "\n";
  Vector3D parent_pos(pos);
  std::array<int, 3> parent_offset(offset);

  for (int i=0; i<3; i++) 
    if (reflection[2-i]) {
      parent_pos[i] *= (-1);
      parent_offset[i] *= (-1);
    }

  return parent->GetNNPotential(parent_pos, parent_offset);
}

double MirrorElectrode::GetNNVoltage(double t, int offset) const {
  return parent->GetNNVoltage(t, offset);
}

void MirrorElectrode::DumpInfo(ostream& os) {
  os << label << "\n";
}

const Vector3D& MirrorElectrode::GetDx() {
  return (*parent).GetDx();
}

void MirrorElectrode::GetFraction(const Vector3D &pos, array<double, 3>& f) {
  (*parent).GetFraction(pos, f);
}
