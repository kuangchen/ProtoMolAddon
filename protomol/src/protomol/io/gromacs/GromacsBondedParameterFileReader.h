/*  -*- c++ -*-  */
#ifndef GROMACSBONDEDPARAMETERFILEREADER_H
#define GROMACSBONDEDPARAMETERFILEREADER_H

#include <protomol/io/Reader.h>
#include <protomol/type/GromacsParameters.h>
 
#include <iostream>
#include <map>

namespace ProtoMol {

  class GromacsBondedParameterFileReader : public Reader {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   public:
    GromacsBondedParameterFileReader();
    explicit GromacsBondedParameterFileReader(const std::string &filename);
    ~GromacsBondedParameterFileReader();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Reader
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   public:
    virtual bool tryFormat();
    virtual bool read();
 
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class GromacsBondedParameterFileReader
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
    bool read(GromacsParameters &gP);

   private:
    void read_bondtypes();
    void read_angletypes();
    void read_dihedrals_and_impropers();
    
    void FillGromacsParameterBondTypes(std::string line);
    void FillGromacsParameterAngleTypes(std::string line);
    void FillGromacsParameterDihedralTypes(std::string line);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   private:
     GromacsParameters *gParams;
 

  };

}

#endif /* GROMACSBONDEDPARAMETERFILEREADER_H */
