#include <protomol/io/PSFReader.h>

#include <protomol/base/StringUtilities.h>
#include <protomol/base/Report.h>
#include <protomol/base/MathUtilities.h>

//____#define DEBUG_PSF

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;
//____PSFReader

PSFReader::PSFReader() :
  Reader(), myPSF(NULL) {}

PSFReader::PSFReader(const string &filename) :
  Reader(filename), myPSF(NULL) {}

PSFReader::~PSFReader() {
  if (myPSF != NULL)
    delete myPSF;
}

bool PSFReader::tryFormat() {
  if (!open())
    return false;
  string psfHead;
  file >> psfHead;
  return file.good() && equalNocase("PSF", psfHead);
}

bool PSFReader::read() {
  if (myPSF == NULL)
    myPSF = new PSF();
  return read(*myPSF);
}

bool PSFReader::read(PSF &psf) {
  if (!tryFormat() || !open()) return false;
  psf.clear();

  // Find header
  string psfHead;
  file >> psfHead;
  if (!equalNocase("PSF", psfHead)) {
    file.setstate(ios::failbit);
    close();
    return false;
  }

  while (!file.eof()) {
    string line(removeBeginEndBlanks(getline()));

    // Exit if nothing more to read
    if (line.empty() && file.eof()) {
      close();
      return true;
    }

    // Find '!', otherwise continue with next line
    if (find(line.begin(), line.end(), '!') == line.end())
      continue;

    // Parse header
    stringstream ss(line);
    vector<string> header;
    int index = -1;
    string str;
    while (ss >> str) {
      if (str[0] == '!' && index < 0)
        index = header.size();
      header.push_back(str);
    }

    if (file.fail() || index < 1 || index > 2 || header.size() < 2) {
      report << recoverable
             << " PSF file \'" << filename
             << "\' has corrupt record header." << endl
             << "\'" << line << "\'"
             << endr;
      file.setstate(ios::failbit);
      return false;
    }
    string keyword = header[index];
    int numrecords = toInt(header[0]);

    // Branch
    if (index == 1) {
      // Normal case: int,keyword

      if (!isInt(header[0])) {
        report << recoverable
               << " PSF file \'" << filename
               <<
        "\' has corrupt record header, number of records should be an integer."
               << endl
               << "\'" << line << "\'"
               << endr;
        file.setstate(ios::failbit);
        close();
        return false;
      }

      if (equalStartNocase("!NTITLE", keyword)) {  // title
        comment = "";
        for (int counter = 0; counter < numrecords; ++counter) {
          // lines of title + one empty line
          line = removeBeginEndBlanks(getline());
          if (!line.empty())
            comment += (comment.empty() ? "" : "\n") + line;
        }

        continue;
      } else if (equalStartNocase("!NATOM", keyword)) {
        for (int counter = 1; counter <= numrecords; ++counter) {
          PSF::Atom temp_atom;
          file >> temp_atom.number  // read atom number
          >> temp_atom.seg_id       // read segment identifier
          >> str                    // read in residue sequence
          >> temp_atom.residue_name // read in residue name
          >> temp_atom.atom_name    // read in atom name
          >> temp_atom.atom_type    // read in name of the second atom
          >> temp_atom.charge       // read in charge
          >> temp_atom.mass         // read in mass
          >> temp_atom.identity;    // atom's chemical identity (used only by
                                    // iSGProtomol)
          temp_atom.residue_sequence = toInt(str);
          if (!isInt(str))
            report << recoverable
                   << "[PSF::read] Expecting a number for residue sequence. "
                   << "I do not know what to do with \'" << str << "\'." <<
            endr;
          psf.atoms.push_back(temp_atom);
          if (file.fail() || (file.eof() && counter < numrecords)) {
            file.setstate(ios::failbit);
            close();
            return false;
          }
        }

        continue;
      } else if (equalStartNocase("!NBOND", keyword)) {
        for (int counter = 1; counter <= numrecords; ++counter) {
          PSF::Bond temp_bond;
          temp_bond.number = counter;
          file >> temp_bond.atom1;
          file >> temp_bond.atom2;
          psf.bonds.push_back(temp_bond);
          if (file.fail() || (file.eof() && counter < numrecords)) {
            file.setstate(ios::failbit);
            close();
            return false;
          }
        }

        continue;
      } else if (equalStartNocase("!NTHETA", keyword)) {
        for (int counter = 1; counter <= numrecords; ++counter) {
          PSF::Angle temp_angle;
          temp_angle.number = counter;
          file >> temp_angle.atom1;
          file >> temp_angle.atom2;
          file >> temp_angle.atom3;
          psf.angles.push_back(temp_angle);
          if (file.fail() || (file.eof() && counter < numrecords)) {
            file.setstate(ios::failbit);
            close();
            return false;
          }
        }

        continue;
      } else if (equalStartNocase("!NPHI", keyword)) {
        for (int counter = 1; counter <= numrecords; ++counter) {
          PSF::Dihedral temp_dihedral;
          temp_dihedral.number = counter;
          file >> temp_dihedral.atom1;
          file >> temp_dihedral.atom2;
          file >> temp_dihedral.atom3;
          file >> temp_dihedral.atom4;
          psf.dihedrals.push_back(temp_dihedral);
          if (file.fail() || (file.eof() && counter < numrecords)) {
            file.setstate(ios::failbit);
            close();
            return false;
          }
        }

        continue;
      } else if (equalStartNocase("!NIMPHI", keyword)) {
        for (int counter = 1; counter <= numrecords; ++counter) {
          PSF::Improper temp_improper;
          temp_improper.number = counter;
          file >> temp_improper.atom1;
          file >> temp_improper.atom2;
          file >> temp_improper.atom3;
          file >> temp_improper.atom4;
          psf.impropers.push_back(temp_improper);
          if (file.fail() || (file.eof() && counter < numrecords)) {
            file.setstate(ios::failbit);
            close();
            return false;
          }
        }

        continue;
      } else if (equalStartNocase("!NDON", keyword)) {
        for (int counter = 1; counter <= numrecords; ++counter) {
          PSF::Donor temp_donor;
          temp_donor.number = counter;
          file >> temp_donor.atom1;
          file >> temp_donor.atom2;
          psf.donors.push_back(temp_donor);
          if (file.fail() || (file.eof() && counter < numrecords)) {
            file.setstate(ios::failbit);
            close();
            return false;
          }
        }

        continue;
      } else if (equalStartNocase("!NACC", keyword)) {
        for (int counter = 1; counter <= numrecords; ++counter) {
          PSF::Acceptor temp_acceptor;
          temp_acceptor.number = counter;
          file >> temp_acceptor.atom1;
          file >> temp_acceptor.atom2;
          psf.acceptors.push_back(temp_acceptor);
          if (file.fail() || (file.eof() && counter < numrecords)) {
            file.setstate(ios::failbit);
            close();
            return false;
          }
        }

        continue;
      } else if (equalStartNocase("!NNB", keyword)) {
        for (int counter = 1; counter <= numrecords; ++counter) {
          PSF::Nonbonded temp_nonbonded;
          temp_nonbonded.number = counter;
          file >> temp_nonbonded.atom1;
          psf.nonbondeds.push_back(temp_nonbonded);
          if (file.fail() || (file.eof() && counter < numrecords)) {
            report << recoverable << "[PSF::read] Expecting " <<
            numrecords << " NNB, read " << counter - 1 <<
            ", reached end of file." << endr;
            close();
            return false;
          }
        }

        continue;
      }
    }

    if (index == 2)
      if (equalStartNocase("!NGRP", keyword)) {
        for (int counter = 1; counter <= numrecords; ++counter) {
          PSF::Ngrp temp_ngrp;
          temp_ngrp.number = counter;
          file >> temp_ngrp.atom1;
          file >> temp_ngrp.atom2;
          file >> temp_ngrp.atom3;
          psf.ngrp.push_back(temp_ngrp);
          if (file.fail() || (file.eof() && counter < numrecords)) {
            report << recoverable << "[PSF::read] Expecting " <<
            numrecords << " NGRP, read " << counter - 1 <<
            ", reached end of file." << endr;
            close();
            return false;
          }
        }

        continue;
      }

    report << recoverable << "[PSF::read] Record " << keyword << " with " <<
    numrecords << " entries not recognized." << endr;
  }
  close();
  return !file.fail();
}

PSF *PSFReader::orphanPSF() {
  PSF *tmp = myPSF;
  myPSF = NULL;
  return tmp;
}

namespace ProtoMol {
  PSFReader &operator>>(PSFReader &psfReader, PSF &psf) {
    psfReader.read(psf);
    return psfReader;
  }
}

