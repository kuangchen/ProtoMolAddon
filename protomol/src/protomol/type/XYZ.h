/*  -*- c++ -*-  */
#ifndef XYZ_H
#define XYZ_H

#include <protomol/type/Vector3DBlock.h>
#include <string>

namespace ProtoMol {
  //_________________________________________________________________XYZ
  /**
   * Container holding coordinates/Vector3D and names
   */
  struct XYZ {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    XYZ() {};
    XYZ(size_t n) :
      coords(n, Vector3D(0.0, 0.0, 0.0)), names(n, std::string()) {}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class XYZ
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    void clear() {coords.clear(); names.clear();}
    size_t size() const {return coords.size();}
    void resize(size_t n) {coords.resize(n), names.resize(n);}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // XYZ container
    Vector3DBlock coords;
    std::vector<std::string> names;
  };
}
#endif /* XYZ_H */
