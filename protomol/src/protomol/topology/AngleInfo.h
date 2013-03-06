/* -*- c++ -*- */
#ifndef ANGLEINFO_H
#define ANGLEINFO_H

#include <vector>
#include <protomol/type/Real.h>

namespace ProtoMol {
  /**
     No comment Cplt!
   */
  class AngleInfo {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Enums & types
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    enum AngleType {ANGLE_VALUE = 10001,
                    ANGLE_POINTER = 10002,
                    ANGLE_NOTSET = 10003,
                    VISITED = 10001,
                    NOT_VISITED = 10002};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    AngleInfo();
    ~AngleInfo();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CubicCellManager
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void addBond(int atom);
    void setAngle(Real angle);
    void setPointer(unsigned int atom);
    void setExclusionAtom();
    void setInnerAtom();
    void setVisited();
    void setAtomID(unsigned int ID);
    AngleType getAngleType() const;
    unsigned int getBond(unsigned int index) const;
    unsigned int numBonds() const;
    unsigned int getAtomID() const;
    Real getAngle() const;
    unsigned int getPointer() const;
    bool isExclusionAtom() const;
    bool isInnerAtom() const;
    bool isVisited() const;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Real my_angle;
    unsigned int my_atomID;
    unsigned int my_pointer;
    AngleType my_angleType;
    bool my_visited;
    bool my_isExclusionAtom;
    bool my_isInnerAtom;
    std::vector<unsigned int> my_bondedAtoms;
  };
}
#endif
