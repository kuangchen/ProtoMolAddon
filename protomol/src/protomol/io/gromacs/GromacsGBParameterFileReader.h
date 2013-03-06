/*  -*- c++ -*-  */
#ifndef GROMACSGBPARAMETERFILEREADER_H
#define GROMACSGBPARAMETERFILEREADER_H


#include <protomol/io/Reader.h>
#include <protomol/type/GromacsParameters.h>
 
#include <iostream>

namespace ProtoMol {

  class GromacsGBParameterFileReader : public Reader {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    public:
      GromacsGBParameterFileReader();
      explicit GromacsGBParameterFileReader(const std::string &filename);
      virtual ~GromacsGBParameterFileReader();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Reader
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    public:
      virtual bool tryFormat();
      virtual bool read();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class GromacsGBParameterFileReader
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    bool read(GromacsParameters &gP);

   private:

    bool read_gromacs_GB_Parameters(std::string line);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   private:
     GromacsParameters *gParams;
       

  };

}

#endif /* GROMACSGBPARAMETERFILEREADER_H */
