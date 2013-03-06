/*  -*- c++ -*-  */
#ifndef GROMACSPARAMETERFILEREADER_H
#define GROMACSPARAMETERFILEREADER_H

#include <protomol/io/Reader.h>
#include <protomol/type/GromacsParameters.h>

#include <iostream>
#include <map>



namespace ProtoMol {

  class GromacsParameterFileReader : public Reader {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   public:
     GromacsParameterFileReader();
     explicit GromacsParameterFileReader(const std::string &filename);
     virtual ~GromacsParameterFileReader();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Reader
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   public:
    virtual bool tryFormat();
    virtual bool read();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class GromacsParameterFileReader
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bool read(GromacsParameters &gP);

   private:
    std::string ParseInclude(std::string line);
    void read_defaults_information();
    void FillGromacsParameterDefaults(std::string s);
  
   public:
    std::string GetBondedParameterFile() { return bondedParameterFile; }
    std::string GetNonbondedParameterFile() { return nonbondedParameterFile; }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   private:
     GromacsParameters *gParams;
     std::string bondedParameterFile;
     std::string nonbondedParameterFile;
   

  };

}

#endif /* GROMACSPARAMETERFILEREADER_H */
