#include <protomol/addon/ToFDetector.h>
#include <protomol/base/PMConstants.h>
#include <protomol/addon/Constants.h>
#include <algorithm>

using namespace ProtoMolAddon;
using namespace ProtoMol;
using namespace ProtoMolAddon::Constant;

ToFDetector::ToFDetector(double detector_pos):
  detector_pos(detector_pos) {

}

ToFDetector::~ToFDetector() {}

void ToFDetector::Initialize(const ProtoMolApp* app) {
  int size = app->positions.size();

  hit.resize(size, false);
  hit_time.resize(size);
  hit_position.resize(size);
  hit_velocity.resize(size);
  atom_name.resize(size);
}

void ToFDetector::Update(const ProtoMolApp* app) {
  for (int i=0; i<app->positions.size(); i++) {
    if (hit[i]) continue;
    else if (app->positions[i][0] * POSITION_CONV > detector_pos)//  &&
      //sqrt(app->positions[i][1] * app->positions[i][1] + app->positions[i][2] * app->positions[i][2]) * POSITION_CONV < 0.8 * 0.0254) {
      {
	hit[i] = true;
	hit_time[i] = app->topology->time * TIME_CONV - 1e-6 + 0.5/1.8e6;
	hit_position[i] = app->positions[i];
	hit_velocity[i] = app->velocities[i];
	atom_name[i] = app->topology->atoms[i].name;
      }
  }
}

namespace ProtoMolAddon {

  ostream& operator<< (ostream& os, ToFDetector& detector) {
    for (int i=0; i<detector.hit_time.size(); i++)
      os << detector.atom_name[i] << "\t" 
	 << detector.hit_time[i] << "\t" 
	 << detector.hit_position[i] * POSITION_CONV << "\t" 
	 << detector.hit_velocity[i] * VELOCITY_CONV << "\n";

    return os;
  }
}
