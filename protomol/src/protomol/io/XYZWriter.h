/*  -*- c++ -*-  */
#ifndef XYZWRITER_H
#define XYZWRITER_H

#include <protomol/io/Writer.h>
#include <protomol/type/XYZ.h>
#include <protomol/topology/Atom.h>
#include <protomol/topology/AtomType.h>

namespace ProtoMol {
  //____XYZWriter

  /**
   * Writes a XYZ format file (ASCII).@n
   *
   * Format:@n
   * 15@n
   * A nice comment!
   * C        1.58890       -1.44870       -0.47000@n
   * C        1.54300       -2.25990        0.77910@n
   * C        2.21440       -3.47410        0.88040@n
   * C        2.16940       -4.22350        2.03080@n
   * ...
   */
  class XYZWriter : public Writer {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    XYZWriter();
    explicit XYZWriter(const std::string &filename);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class XYZ
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    bool write(const XYZ &xyz);
    bool write(const Vector3DBlock &coords,
               const std::vector<std::string> &names);
    bool write(const Vector3DBlock &coords, const std::vector<Atom> &atoms,
               const std::vector<AtomType> &atomTypes);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Friends
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    friend XYZWriter &operator<<(XYZWriter &xyzWriter, const XYZ &xyz);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
  };

  //____INLINES
}
#endif /* XYZWRITER_H */
