/*  -*- c++ -*-  */
#ifndef PDBWRITER_H
#define PDBWRITER_H

#include <protomol/io/Writer.h>
#include <protomol/type/PDB.h>

namespace ProtoMol {
  //____PDBWriter

  /**
   * Writes a PDB (ASCII) file, ATOM record only.
   */
  class PDBWriter : public Writer {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    PDBWriter();
    explicit PDBWriter(const std::string &filename);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class PDB
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    bool write(const PDB &pdb);
    bool write(const Vector3DBlock &coords, const PDB &pdb);
    bool write(const Vector3DBlock &coords,
               const std::vector<PDB::Atom> &atoms);
    bool write(const Vector3DBlock &coords, const std::vector<PDB::Atom> &atoms,
               const std::vector<PDB::Ter> &ters);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Friends
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    friend PDBWriter &operator<<(PDBWriter &pdbWriter, const PDB &pdb);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
  };

  //____INLINES
}
#endif /* PDBWRITER_H */
