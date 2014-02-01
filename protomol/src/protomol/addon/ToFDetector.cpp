#include <protomol/addon/ToFDetector.h>
#include <protomol/base/PMConstants.h>
#include <protomol/addon/Constants.h>
#include <algorithm>

using namespace ProtoMolAddon;
using namespace ProtoMol;
using namespace ProtoMolAddon::Constant;

ToFDetector::AtomRecord::AtomRecord(const string &name, 
				    const Vector3D &pos, 
				    const Vector3D &vel, 
				    double time) : 
  name(name), pos(pos), vel(vel), time(time), hit(false) {
}

ostream& ProtoMolAddon::operator<< (ostream &os, const ToFDetector::AtomRecord &r) {
  os << r.name << "\t" 
     << r.time << "\t" 
     << r.pos * POSITION_CONV << "\t" 
     << r.vel * VELOCITY_CONV << std::endl;

  return os;
}


ToFDetector::ToFDetector(double detector_pos):
  record_list(), detector_pos(detector_pos) 
{}


void ToFDetector::Initialize(const ProtoMolApp *app) {
  for (size_t i=0; i<app->positions.size(); i++) 
    record_list.push_back(shared_ptr<AtomRecord>(new AtomRecord(app->topology->atoms[i].name,
								app->positions[i],
								app->velocities[i])));
}


void ToFDetector::UpdateRecord(const ProtoMolApp *app) {
  for (auto &r: record_list) {
    if (r->Hit()) continue;
    else if (r->pos[0] * POSITION_CONV > detector_pos)
      r->SetHitTime(app->topology->time * TIME_CONV);
  }
}


ostream& ProtoMolAddon::operator<< (ostream &os, const ToFDetector &detector) {
  for (auto &r : detector.record_list) os << r;

  return os;
}

