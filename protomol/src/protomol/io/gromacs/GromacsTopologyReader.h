/*  -*- c++ -*-  */
#ifndef GROMACSTOPOLOGYREADER_H
#define GROMACSTOPOLOGYREADER_H


#include <protomol/io/Reader.h>
#include <protomol/type/GromacsTopology.h>

#include <iostream>
#include <map>

namespace ProtoMol {

  //___________________________________________ GromacsTopologyReader

  class GromacsTopologyReader : public Reader {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   public:
    GromacsTopologyReader();
    explicit GromacsTopologyReader(const std::string &filename);
    virtual ~GromacsTopologyReader();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Reader
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   public:
    virtual bool tryFormat();
    virtual bool read();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class GromacsTopologyReader
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bool read(GromacsTopology &gT);
  
    std::string GetParameterFilename();


   private:
    std::string ParseInclude(std::string s);
    void read_atom_information();
    void read_bond_information();
    void read_angle_information();
    void read_dihedral_information();
    void FillGromacsTopologyAtom(std::string line);
    void FillGromacsTopologyBond(std::string line);
    void FillGromacsTopologyAngle(std::string line);
    void FillGromacsTopologyDihedrals(std::string line);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   private:

    GromacsTopology *gTopo;

    std::vector<std::string> parameter_files;

  };

}

#endif /* GROMACSTOPOLOGYREADER_H */
