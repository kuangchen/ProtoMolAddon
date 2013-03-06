/*  -*- c++ -*-  */
#ifndef XYZTRAJECTORYREADER_H
#define XYZTRAJECTORYREADER_H

#include <protomol/io/XYZReader.h>
#include <protomol/type/XYZ.h>

namespace ProtoMol {
  //____XYZTrajectoryReader

  /**
   * Reads a XYY trajectory files (ASCII).
   */
  class XYZTrajectoryReader : public XYZReader {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    XYZTrajectoryReader();
    explicit XYZTrajectoryReader(const std::string &filename);
    virtual ~XYZTrajectoryReader();

    virtual bool tryFormat();
    virtual bool read();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class XYZTrajectoryReader
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bool read(Vector3DBlock &xyz);
    void doRead(Vector3DBlock &xyz);

    Vector3DBlock *orphanXYZ();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Friends
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    friend XYZTrajectoryReader &operator>>(XYZTrajectoryReader &reader,
                                           Vector3DBlock &xyz);

  private:
    Vector3DBlock* xyz;
    bool first;
    //std::vector<XYZ> *xyz;
  };
}
#endif /* XYZTRAJECTORYREADER_H */
