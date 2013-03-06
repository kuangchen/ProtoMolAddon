/*  -*- c++ -*-  */
#ifndef PDB_H
#define PDB_H

#include <string>
#include <vector>

#include <protomol/type/Vector3DBlock.h>
#include <protomol/base/Report.h>

namespace ProtoMol {
  //___________________________________________________________________PDB
  /**
   * PDB container holding Atom's and coordinates
   */
  struct PDB {
    //___________________________________________________________________Atom
    /**
     * This class holds data for a basic PDB element; excluding the
     * coordinates (x,y,z), stored separately.
     */
    struct Atom {
      Atom();
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // My data members
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Atom(std::string elementType,
           int elementNum,
           std::string elementName,
           std::string altLoc,
           std::string residueName,
           std::string chainID,
           int residueNum,
           std::string insertionCode,
           Real occupancy,
           Real tempFactor,
           std::string segID,
           std::string symbol,
           std::string charge,
           int hvyAtomGrpsize);

      std::string elementType;     ///< record_name
      int elementNum;      ///< serial_size
      std::string elementName;     ///< atom_name
      std::string altLoc;      ///< alternate_location
      std::string residueName;     ///< residue_name
      std::string chainID;     ///< chain_id
      int residueNum;      ///< residue_sequence
      std::string insertionCode;   ///< insertion_code
      //Real x;
      //Real y;
      //Real z;
      Real occupancy;       ///< occupancy
      Real tempFactor;      ///< temp_factor
      std::string segID;       ///< seg_id
      std::string symbol;      ///< element_symbol
      std::string charge;      ///< charge
      int hvyAtomGrpsize;  ///< ??? throw in zeros

      /// Record element start positions
      enum Start {
        S_RECORD_NAME = 1 - 1,
        S_SERIAL = 7 - 1,
        S_ATOM_NAME = 13 - 1,
        S_ALT_LOC = 17 - 1,
        S_RES_NAME = 18 - 1,
        S_CHAIN_ID = 22 - 1,
        S_RES_SEQ = 23 - 1,
        S_I_CODE = 27 - 1,
        S_X = 31 - 1,
        S_Y = 39 - 1,
        S_Z = 47 - 1,
        S_OCCUP = 55 - 1,
        S_TEMP_FACT = 61 - 1,
        S_SEG_ID = 73 - 1,
        S_ELEMENT_SYMBOL = 77 - 1,
        S_CHARGE = 79 - 1
      };

      /// Record element length
      enum Length {
        L_RECORD_NAME = 1 - (1 - 6),
        L_SERIAL = 1 - (7 - 11),
        L_ATOM_NAME = 1 - (13 - 16),
        L_ALT_LOC = 1 - (17 - 17),
        L_RES_NAME = 1 - (18 - 20),
        L_CHAIN_ID = 1 - (22 - 22),
        L_RES_SEQ = 1 - (23 - 26),
        L_I_CODE = 1 - (27 - 27),
        L_X = 1 - (31 - 38),
        L_Y = 1 - (39 - 46),
        L_Z = 1 - (47 - 54),
        L_OCCUP = 1 - (55 - 60),
        L_TEMP_FACT = 1 - (61 - 66),
        L_SEG_ID = 1 - (73 - 76),
        L_ELEMENT_SYMBOL = 1 - (77 - 78),
        L_CHARGE = 1 - (79 - 80)
      };


      friend Report::MyStreamer &operator<<(Report::MyStreamer &OS,
                                            const Atom &p);
    };

    struct Ter {
      Ter();
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // My data members
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Ter(std::string elementType,
          int elementNum,
          std::string residueName,
          std::string chainID,
          int residueNum,
          std::string insertionCode);

      std::string elementType;     ///< record_name
      int elementNum;      ///< serial_size
      std::string residueName;     ///< residue_name
      std::string chainID;     ///< chain_id
      int residueNum;      ///< residue_sequence
      std::string insertionCode;   ///< insertion_code

      /// Record element start positions
      enum Start {
        S_RECORD_NAME = 1 - 1,
        S_SERIAL = 7 - 1,
        S_RES_NAME = 18 - 1,
        S_CHAIN_ID = 22 - 1,
        S_RES_SEQ = 23 - 1,
        S_I_CODE = 27 - 1
      };

      /// Record element length
      enum Length {
        L_RECORD_NAME = 1 - (1 - 6),
        L_SERIAL = 1 - (7 - 11),
        L_RES_NAME = 1 - (18 - 20),
        L_CHAIN_ID = 1 - (22 - 22),
        L_RES_SEQ = 1 - (23 - 26),
        L_I_CODE = 1 - (27 - 27)
      };
    };

    //___________________________________________________________________PDB

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PDB() {};
    PDB(size_t n) : coords(n, Vector3D(0.0, 0.0, 0.0)), atoms(n) {}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class PDB
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    void clear();
    size_t size() const {return coords.size();}
    void resize(size_t n) {coords.resize(n), atoms.resize(n);}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Vector3DBlock coords;
    std::vector<PDB::Atom> atoms;
    std::vector<PDB::Ter> ters;
  };


  //___________________________________________________________________INLINES
}

#endif /* PDB_H */




