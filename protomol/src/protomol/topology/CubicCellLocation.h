/* -*- c++ -*- */
#ifndef CUBICCELLLOCATION_H
#define CUBICCELLLOCATION_H

namespace ProtoMol {
  //________________________________________ CubicCellLocation
  /**
   * Equal-sized (cubic) cell.
   */
  class CubicCellLocation {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    CubicCellLocation(void) : x(0), y(0), z(0) {}
    CubicCellLocation(int a, int b, int c) : x(a), y(b), z(c) {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CubicCellLocation
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    bool operator<(const CubicCellLocation &) const;
    bool operator==(const CubicCellLocation &c) const {
      return x == c.x && y == c.y && z == c.z;
    }
    CubicCellLocation operator-(const CubicCellLocation &b) const {
      return CubicCellLocation(x - b.x, y - b.y, z - b.z);
    }
    CubicCellLocation operator+(const CubicCellLocation &b) const {
      return CubicCellLocation(x + b.x, y + b.y, z + b.z);
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    int x;
    int y;
    int z;
  };

  //________________________________________ INLINES

  inline bool CubicCellLocation::operator<(const CubicCellLocation &c) const {
    if (x < c.x) return true;
    if (x > c.x) return false;
    if (y < c.y) return true;
    if (y > c.y) return false;
    if (z < c.z) return true;
    return false;
  }
}
#endif /* CUBICCELLLOCATION_H */
