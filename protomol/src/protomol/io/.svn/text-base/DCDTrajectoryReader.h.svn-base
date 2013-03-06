/*  -*- c++ -*-  */
#ifndef DCDTRAJECTORYREADER_H
#define DCDTRAJECTORYREADER_H

#include <protomol/io/Reader.h>
#include <protomol/type/XYZ.h>
#include <protomol/type/TypeSelection.h>

#include <vector>
#include <cstring>

namespace ProtoMol {
  //____DCDTrajectoryReader

  /**
   * Reads a DCD trajectory file, frame by frame. Automatic endianess
   * detection.
   */
  class DCDTrajectoryReader : public Reader {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Typedef
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    typedef TypeSelection::Int<4>::type int32;
    typedef TypeSelection::Float<4>::type float4;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    DCDTrajectoryReader();
    explicit DCDTrajectoryReader(const std::string &filename);
    virtual ~DCDTrajectoryReader();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Reader
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool tryFormat();
    virtual bool read();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class DCDTrajectory
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    bool read(Vector3DBlock &coords);
    void doRead(Vector3DBlock &coords);
    int readFirstStep();

    Vector3DBlock *orphanXYZ();

  private:
    void fortranRead(char *data, unsigned int size,
                     const std::string &err = "");
    char *fortranReadX(char *data, unsigned int &size,
                      const std::string &err = "");

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Friends
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    friend DCDTrajectoryReader &operator>>(DCDTrajectoryReader &reader,
                                           Vector3DBlock &xyz);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Vector3DBlock* xyz;
    //std::vector<XYZ> *xyz;
  private:
    bool swap;
    bool first;
    int32 natoms;
    // Read header
    struct {
      char cord[4];
      int32 frames;
      int32 firststep;
      char ignore1[24];
      int32 freeIndexes;
      char ignore2[44];
    } header;
  };
}
#endif /* DCDTRAJECTORYREADER_H */

