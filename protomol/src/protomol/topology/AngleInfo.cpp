#include <protomol/topology/AngleInfo.h>

using namespace ProtoMol;
AngleInfo::AngleInfo() {
  my_angle = 0.0;
  my_pointer = 0;
  my_angleType = ANGLE_NOTSET;
  my_visited = false;
  my_isExclusionAtom = false;
  my_isInnerAtom = false;
  my_atomID = 0;
}

AngleInfo::~AngleInfo() {}

void AngleInfo::setAtomID(unsigned int ID) {
  my_atomID = ID;
}

unsigned int AngleInfo::getAtomID() const {
  return my_atomID;
}

void AngleInfo::addBond(int atom) {
  my_bondedAtoms.push_back(atom);
}

void AngleInfo::setAngle(Real angle) {
  my_angle = angle;
  my_angleType = ANGLE_VALUE;
}

void AngleInfo::setPointer(unsigned int atom) {
  my_pointer = atom;
  my_angleType = ANGLE_POINTER;
}

AngleInfo::AngleType AngleInfo::getAngleType() const {
  return my_angleType;
}

Real AngleInfo::getAngle() const {
  return my_angle;
}

unsigned int AngleInfo::getPointer() const {
  return my_pointer;
}

void AngleInfo::setVisited() {
  my_visited = true;
}

bool AngleInfo::isVisited() const {
  return my_visited;
}

void AngleInfo::setInnerAtom() {
  my_isInnerAtom = true;
}

bool AngleInfo::isInnerAtom() const {
  return my_isInnerAtom;
}

void AngleInfo::setExclusionAtom() {
  my_isExclusionAtom = true;
}

bool AngleInfo::isExclusionAtom()  const {
  return my_isExclusionAtom;
}

unsigned int AngleInfo::numBonds() const {
  return my_bondedAtoms.size();
}

unsigned int AngleInfo::getBond(unsigned int index) const {
  if (index < my_bondedAtoms.size())
    return my_bondedAtoms[index];
  return 0;
}

