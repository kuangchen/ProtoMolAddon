/*  -*- c++ -*-  */
#ifndef EIGENVECTOTTEXTREADER_H
#define EIGENVECTOTTEXTREADER_H

#include <protomol/io/Reader.h>
#include <protomol/type/EigenvectorInfo.h>
#include <protomol/type/TypeSelection.h>

namespace ProtoMol {
  //____PDBReader

  /**
   * Reads an EigenvectorInfo ASCII file.
   */
  class EigenvectorTextReader : public Reader {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    EigenvectorTextReader();
    explicit EigenvectorTextReader(const std::string &filename);
    explicit EigenvectorTextReader(const char *filename);
    virtual ~EigenvectorTextReader();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class File
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool open() {return File::open();}

    virtual bool open(const std::string &filename) {
      return File::open(filename);
    }

    virtual bool open(const char *filename) {return File::open(filename);}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Reader
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool tryFormat() {return true;}

    virtual bool read();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class Eigenvector
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    bool read(EigenvectorInfo &ei);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Friends
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    friend EigenvectorTextReader &operator>>
    (EigenvectorTextReader &eigenvectorReader, EigenvectorInfo &info);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    typedef TypeSelection::Int<4>::type int32;
    EigenvectorInfo *myEigenvectorInfo;
  };
}
#endif /* EIGENVECTOTTEXTREADER_H */
