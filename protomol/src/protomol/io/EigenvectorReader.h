/*  -*- c++ -*-  */
#ifndef EIGENVECTORREADER_H
#define EIGENVECTORREADER_H

#include <protomol/io/Reader.h>
#include <protomol/type/EigenvectorInfo.h>
#include <protomol/type/TypeSelection.h>

namespace ProtoMol {
  //_________________________________________________________________PDBReader
  /**
   * Reads an EigenvectorInfo binary file.
   */
  class EigenvectorReader : public Reader {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    EigenvectorReader();
    explicit EigenvectorReader(const std::string &filename);
    explicit EigenvectorReader(const char *filename);
    virtual ~EigenvectorReader();

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
    virtual bool read(EigenvectorInfo &ei);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Friends
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    friend EigenvectorReader &operator>>(EigenvectorReader &eigenvectorReader,
                                         EigenvectorInfo &info);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    typedef TypeSelection::Int<4>::type int32;
    EigenvectorInfo *myEigenvectorInfo;
    bool mySwapEndian;
  };
}
#endif /* EIGENVECTORREADER_H */
