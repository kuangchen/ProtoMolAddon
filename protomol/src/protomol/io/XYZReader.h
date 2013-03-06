/*  -*- c++ -*-  */
#ifndef XYZREADER_H
#define XYZREADER_H

#include <protomol/io/Reader.h>
#include <protomol/type/XYZ.h>

#ifdef HAVE_LIBFAH
namespace FAH {
  class Exception;
};
typedef FAH::Exception Exception;
#endif

namespace ProtoMol {
  //_________________________________________________________________XYZReader
  /**
   * Reads a XYZ format file (ASCII).
   */
  class XYZReader : public Reader {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    XYZReader();
    explicit XYZReader(const std::string &filename);
    virtual ~XYZReader();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Reader
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool tryFormat();
    virtual bool read();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class XYZ
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    bool read(XYZ &xyz);
    bool read(Vector3DBlock &coords, std::vector<std::string> &names);
    void doRead(Vector3DBlock &coords, std::vector<std::string> &names);

    XYZ getXYZ() const;
    Vector3DBlock *orphanCoords();
    std::vector<std::string> *orphanNames();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Friends
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    friend XYZReader &operator>>(XYZReader &xyzReader, XYZ &xyz);
    friend XYZReader &operator>>(XYZReader &xyzReader, Vector3DBlock &coords);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Vector3DBlock *myCoords;
    std::vector<std::string> *myNames;
  };
}
#endif /* XYZREADER_H */
