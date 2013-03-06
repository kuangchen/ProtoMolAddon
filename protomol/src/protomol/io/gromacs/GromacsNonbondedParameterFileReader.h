/*  -*- c++ -*-  */
#ifndef GROMACSNONBONDEDPARAMETERFILEREADER_H
#define GROMACSNONBONDEDPARAMETERFILEREADER_H

#include <protomol/io/Reader.h>
#include <protomol/type/GromacsParameters.h>

#include <iostream>
#include <map>

namespace ProtoMol {

  class GromacsNonbondedParameterFileReader : public Reader {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   public:
    GromacsNonbondedParameterFileReader();
    explicit GromacsNonbondedParameterFileReader(const std::string &filename);
    ~GromacsNonbondedParameterFileReader();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Reader
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   public:
    virtual bool tryFormat();
    virtual bool read();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class GromacsNonbondedParameterFileReader
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    bool read(GromacsParameters &gP);

   private:
    void read_atomtypes();
    void FillGromacsParameterAtomTypes(std::string s);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   private:

    GromacsParameters *gParams;

  };

}

#endif /* GROMACSNONBONDEDPARAMETERFILEREADER_H */
