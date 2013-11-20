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
}

void ToFDetector::Update(const ProtoMolApp* app) {
  for (int i=0; i<app->positions.size(); i++) {
    if (hit[i]) continue;
    else if (app->positions[i][0] * POSITION_CONV > detector_pos) {
      hit[i] = true;
      hit_time[i] = app->topology->time * TIME_CONV - 1e-6;
    }
  }
}

