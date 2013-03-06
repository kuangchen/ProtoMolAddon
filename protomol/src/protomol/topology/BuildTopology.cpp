#include <protomol/topology/BuildTopology.h>

#include <protomol/base/Exception.h>
#include <protomol/type/PAR.h>
#include <protomol/type/PSF.h>
#include <protomol/type/String.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/CoulombSCPISMParameterTable.h>
#include <protomol/topology/LennardJonesParameterTable.h>

#include <vector>

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;


//____findNextNeighbor
static void findNextNeighbor(int a, vector<int> &v, vector<PairInt> &p,
                             vector<char> &unused,
                             const vector<vector<int> > &graph,
                             set<PairInt> &pairs) {
  if (unused[a] > 0) {
    v.push_back(a);
    unused[a] = 0;
    for (unsigned int i = 0; i < graph[a].size(); i++) {
      set<PairInt>::iterator itr =
        pairs.find(PairInt(min(a, graph[a][i]), max(a, graph[a][i])));

      if (itr != pairs.end()) {
        p.push_back(PairInt(a, graph[a][i]));
        pairs.erase(itr);
      }
      findNextNeighbor(graph[a][i], v, p, unused, graph, pairs);
    }
  }
}

void ProtoMol::buildTopology(GenericTopology *topo, const PSF &psf,
                             const PAR &par, bool dihedralMultPSF,  
                             CoulombSCPISMParameterTable *mySCPISMTable) { 
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // First, generate the array of atomtypes
  // Each time a new atom comes up, we need to check if it is
  // already in the vector....
  // NOTE:  this may take a while for large systems; however, it will cut
  // down on the size of the atomTypes vector, and therefore, the amount
  // access time in the back end.
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  topo->atoms.clear();
  topo->atomTypes.clear();
  topo->bonds.clear();
  topo->angles.clear();
  topo->dihedrals.clear();
  topo->impropers.clear();

  //Ryckert-Belleman
  topo->rb_dihedrals.clear();

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Get the atoms
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  map<string, int> atomLookUpTable;

  // loop over all atoms in the PSF object
  for (vector<PSF::Atom>::const_iterator atom = psf.atoms.begin();
       atom != psf.atoms.end(); ++atom) {
    // Two data members for AtomType, name and mass
    Atom tempatom;
    AtomType tempatomtype;
    tempatomtype.name = atom->atom_type;
    tempatomtype.mass = atom->mass;
    tempatomtype.symbolName = atomTypeToSymbolName(atom->atom_type);
    tempatomtype.charge = atom->charge;

    // Now check if this already exists (same name)
    if (atomLookUpTable.find(tempatomtype.name) == atomLookUpTable.end()) {
      atomLookUpTable[tempatomtype.name] = topo->atomTypes.size();
      topo->atomTypes.push_back(tempatomtype);
    }

    // First, we need to find the index. (an integer corresponding
    // to the type of the atom
    tempatom.name = atom->atom_name;
    tempatom.type = atomLookUpTable[tempatomtype.name];
    tempatom.residue_name = atom->residue_name;
    tempatom.residue_seq = atom->residue_sequence;
    // Now, the scaled charge.  This is straightforward.
    tempatom.scaledCharge = (atom->charge) * Constant::SQRTCOULOMBCONSTANT;
    tempatom.scaledMass = atom->mass;
    // Now we need the size of the group for heavy atom ordering
    // We need to parse the name for any H's then any numbers following
    // First, if the atom is an H then this is 0
    if (atom->atom_type == "H") tempatom.hvyAtom = 0;
    else {
      // Otherwise, we need to parse..
      // Initialize to 1
      tempatom.hvyAtom = 1;
      for (unsigned int pos = 0; pos < atom->atom_type.size(); ++pos)
        if (atom->atom_type[pos] == 'H') {
          string number = "";
          while (isdigit(atom->atom_type[++pos]))
            number += atom->atom_type[pos];

          if (number == "")  // never entered loop, default is 1
            number = "1";
          tempatom.hvyAtom += atoi(number.c_str());
        }
    }
    // C/C++ starts at 0, where PSF/PDB at 1
    tempatom.atomNum = atom->number - 1;
    // Also the molecule - using residue sequence for now
    topo->atoms.push_back(tempatom);
  }

  // calculate the # of degrees of freedom, if there are any bond constraints
  // they will be subtracted later by ModifierShake
  topo->degreesOfFreedom = 3 * topo->atoms.size() - 3;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Get the bonds
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // First create look-up-table
  map<string, vector<PAR::Bond>::const_iterator> bondLookUpTable;
  for (vector<PAR::Bond>::const_iterator bond = par.bonds.begin();
       bond != par.bonds.end();
       ++bond)
    bondLookUpTable[String::toUpper(bond->atom1) + "," + String::toUpper(bond->atom2)] = bond;

  // Find the parameters from PAR
  int ignoredBonds = 0;
  for (vector<PSF::Bond>::const_iterator bond = psf.bonds.begin();
       bond != psf.bonds.end(); ++bond) {
    // store the ID numbers of the bonded atoms
    int atom1 = bond->atom1 - 1;
    int atom2 = bond->atom2 - 1;

    // store the type names of the bonded atoms
    string bond1(topo->atomTypes[topo->atoms[atom1].type].name);
    string bond2(topo->atomTypes[topo->atoms[atom2].type].name);

    map<string,
        vector<PAR::Bond>::const_iterator>::const_iterator currentbond =
      bondLookUpTable.find(String::toUpper(bond1) + "," + String::toUpper(bond2));

    if (currentbond == bondLookUpTable.end())
      currentbond = bondLookUpTable.find(String::toUpper(bond2) + "," + String::toUpper(bond1));

    // if we still have not found this bond type in the PAR object, report an
    // error
    if (currentbond == bondLookUpTable.end()) {
      ostringstream err;

      err << "Could not find bond '" << bond1 << "'-'" << bond2 << "' ("
          << bond->atom1 << "," << bond->atom2 << ")" << endl;

      for (map<string, vector<PAR::Bond>::const_iterator>::const_iterator i =
             bondLookUpTable.begin(); i != bondLookUpTable.end(); i++)
        err << i->first << endl;

      THROW(err.str());
    }

    // if we have found this bond type then copy the bond parameters
    // into the topology
    Bond tempbond;
    tempbond.springConstant = currentbond->second->forceConstant;
    tempbond.restLength = currentbond->second->distance;
    tempbond.atom1 = atom1;
    tempbond.atom2 = atom2;
    topo->bonds.push_back(tempbond);

    // populate the vector of bonds maintained at each atom
    topo->atoms[atom1].mybonds.push_back((topo->bonds.size()) - 1);
    topo->atoms[atom2].mybonds.push_back((topo->bonds.size()) - 1);

    if (tempbond.springConstant == 0.0) ignoredBonds++;
  }

  if (ignoredBonds > 0)
    report << hint << "Systems contains " << ignoredBonds
           << " bonds with zero force constants." << endr;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Get the angles
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // First create look-up-table
  map<string, vector<PAR::Angle>::const_iterator> angleLookUpTable;
  for (vector<PAR::Angle>::const_iterator angle = par.angles.begin();
       angle != par.angles.end(); ++angle)
    angleLookUpTable[String::toUpper(angle->atom1) + "," + String::toUpper(angle->atom2) + "," +
                     String::toUpper(angle->atom3)] = angle;

  // Find the parameters from PAR
  int ignoredAngles = 0;

  // loop over the angle list in the PSF object
  for (vector<PSF::Angle>::const_iterator angle = psf.angles.begin();
       angle != psf.angles.end(); ++angle) {
    // store the ID numbers of the atoms in this angle
    int atom1 = angle->atom1 - 1;
    int atom2 = angle->atom2 - 1;
    int atom3 = angle->atom3 - 1;

    // store the type names of the atoms in this angle
    string angle1(topo->atomTypes[topo->atoms[atom1].type].name);
    string angle2(topo->atomTypes[topo->atoms[atom2].type].name);
    string angle3(topo->atomTypes[topo->atoms[atom3].type].name);

    map<string, vector<PAR::Angle>::const_iterator>::const_iterator
    currentangle =
      angleLookUpTable.find(angle1 + "," + angle2 + "," + angle3);

    if (currentangle == angleLookUpTable.end())
      currentangle = angleLookUpTable.find(
        angle3 + "," + angle2 + "," + angle1);

    // if we still have not found this angle type in the PAR object, report an
    // error
    if (currentangle == angleLookUpTable.end())
      THROWS("Could not find angle '" << angle1 << "'-'"
                                      << angle2 << "'-'" << angle3 << "'.");

    // if we have found this angle type then copy the angle parameters
    // into the topology
    Angle tempangle;
    tempangle.atom1 = atom1;
    tempangle.atom2 = atom2;
    tempangle.atom3 = atom3;
    tempangle.forceConstant = currentangle->second->forceConstant;
    tempangle.restAngle = dtor(currentangle->second->angleval);
    if (currentangle->second->ub_flag) {       // do we want defaults for these
      tempangle.ureyBradleyConstant = currentangle->second->k_ub;
      tempangle.ureyBradleyRestLength = currentangle->second->r_ub;
    } else { // no Urey-Bradley term specified
      tempangle.ureyBradleyConstant = 0.0;
      tempangle.ureyBradleyRestLength = 0.0;
    }
    topo->angles.push_back(tempangle);
    if (tempangle.forceConstant == 0.0)
      ignoredAngles++;
  }

  if (ignoredAngles > 0)
    report << hint << "Systems contains " << ignoredAngles <<
    " angles with zero force constants." << endr;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Get the dihedrals
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  // One change I made was to assume that a dihedral will only appear
  // once in the .psf file regardless of it's multiplicity.  The
  // multiplicity should be handled in the .par file.
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // First create look-up-table
  map<string, vector<PAR::Dihedral>::const_iterator> dihedralLookUpTable;
  for (vector<PAR::Dihedral>::const_iterator dihedral = par.dihedrals.begin();
       dihedral != par.dihedrals.end(); ++dihedral)
    dihedralLookUpTable[String::toUpper(dihedral->atom1) + "," + String::toUpper(dihedral->atom2) + "," +
                        String::toUpper(dihedral->atom3) + "," + String::toUpper(dihedral->atom4)] = dihedral;

  // Find the parameters from PAR
  // loop over the dihedral list in the PSF object
  for (vector<PSF::Dihedral>::const_iterator dihedral = psf.dihedrals.begin();
       dihedral != psf.dihedrals.end(); ++dihedral) {
    // store the ID numbers of the atoms in this dihedral
    int atom1 = dihedral->atom1 - 1;
    int atom2 = dihedral->atom2 - 1;
    int atom3 = dihedral->atom3 - 1;
    int atom4 = dihedral->atom4 - 1;

    // store the type names of the atoms in this dihedral
    string dihedral1 = topo->atomTypes[topo->atoms[atom1].type].name;
    string dihedral2 = topo->atomTypes[topo->atoms[atom2].type].name;
    string dihedral3 = topo->atomTypes[topo->atoms[atom3].type].name;
    string dihedral4 = topo->atomTypes[topo->atoms[atom4].type].name;

    map<string,
        vector<PAR::Dihedral>::const_iterator>::const_iterator
    currentdihedral =
      dihedralLookUpTable.find(
        dihedral1 + "," + dihedral2 + "," + dihedral3 + "," + dihedral4);

    // if this dihedral type has not been found, try reversing the order of
    // the atom types
    if (currentdihedral == dihedralLookUpTable.end())
      currentdihedral = dihedralLookUpTable.find(
        dihedral4 + "," + dihedral3 + "," + dihedral2 + "," + dihedral1);

    // Try wildcards if necessary
    if (currentdihedral == dihedralLookUpTable.end()) {
      currentdihedral = dihedralLookUpTable.find(
        "X," + dihedral2 + "," + dihedral3 + ",X");
      if (currentdihedral == dihedralLookUpTable.end())
        currentdihedral = dihedralLookUpTable.find(
          "X," + dihedral3 + "," + dihedral2 + ",X");

      //for GROMACS
      if (currentdihedral == dihedralLookUpTable.end())
         currentdihedral = dihedralLookUpTable.find(
         "X," + dihedral2 + "," + dihedral3 + "," + dihedral4); 

      if (currentdihedral == dihedralLookUpTable.end())
         currentdihedral = dihedralLookUpTable.find(
         "X," + string("X,") + dihedral3 + "," + dihedral4);
    }

    // if we still have not found this dihedral type in the PAR object, report
    // an error
    if (currentdihedral == dihedralLookUpTable.end())
      THROWS(
        "Could not find dihedral '" << dihedral1 << "'-'" << dihedral2
                                    << "'-'" << dihedral3 << "'-'" <<
        dihedral4 << "'.");

    // if we have found this dihedral type then copy the
    // dihedral parameters into the topology
    Torsion torsion;
    torsion.atom1 = atom1;
    torsion.atom2 = atom2;
    torsion.atom3 = atom3;
    torsion.atom4 = atom4;

    torsion.periodicity = currentdihedral->second->periodicity;
    torsion.forceConstant = currentdihedral->second->forceConstant;
    torsion.phaseShift = dtor(currentdihedral->second->phaseShift);
    torsion.multiplicity = currentdihedral->second->multiplicity;
    if (topo->dihedrals.empty() ||
        topo->dihedrals[topo->dihedrals.size() - 1].atom1 != atom1 ||
        topo->dihedrals[topo->dihedrals.size() - 1].atom2 != atom2 ||
        topo->dihedrals[topo->dihedrals.size() - 1].atom3 != atom3 ||
        topo->dihedrals[topo->dihedrals.size() - 1].atom4 != atom4) {
      if (dihedralMultPSF) {
        torsion.periodicity.resize(1);
        torsion.forceConstant.resize(1);
        torsion.phaseShift.resize(1);
        torsion.multiplicity = 1;
      }
      topo->dihedrals.push_back(torsion);
    } else if (dihedralMultPSF) {
      Torsion &tmp = topo->dihedrals[topo->dihedrals.size() - 1];
      if (tmp.multiplicity > torsion.multiplicity)
        THROWS("PSF multiplicity definition of dihedral (" << dihedral1 << ","
               << dihedral2 << "," << dihedral3 << "," << dihedral4 
               << ") exceeded PAR definition.");

        tmp.periodicity.push_back(torsion.periodicity[tmp.multiplicity]);
      tmp.forceConstant.push_back(torsion.forceConstant[tmp.multiplicity]);
      tmp.phaseShift.push_back(torsion.phaseShift[tmp.multiplicity]);
      tmp.multiplicity++;
    } else
      THROWS("Unexpected PSF multiplicity definition of dihedral ("
             << dihedral1 << "," << dihedral2 << "," << dihedral3 << ","
             << dihedral4 << ") occurred, use dihedral multiplicity PSF.");
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Get the impropers
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  // One change I made was to assume that a improper will only appear
  // once in the .psf file regardless of it's multiplicity.  The
  // multiplicity should be handled in the .par file.
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //         No wildcard usage is allowed for bonds and angles. For dihedrals,
  // two types are allowed; A - B - C - D (all four atoms specified) and
  // X - A - B - X (only middle two atoms specified). Double dihedral
  // specifications may be specified for the four atom type by listing a
  // given set twice. When specifying this type in the topology file, specify
  // a dihedral twice (with nothing intervening) and both forms will be used.
  //
  //         There are five choices for wildcard usage for improper dihedrals;
  // 1) A - B - C - D  (all four atoms, double specification allowed)
  // 2) A - X - X - B
  // 3) X - A - B - C
  // 4) X - A - B - X
  // 5) X - X - A - B
  // When classifying an improper dihedral, the first acceptable match (from
  // the above order) is chosen. The match may be made in either direction
  // ( A - B - C - D = D - C - B - A).
  //
  //         The periodicity value for dihedrals and improper dihedral terms
  // must be an integer. If it is positive, then a cosine functional form is
  // used. Only positive values of 1,2,3,4,5 and 6 are allowed for the vector,
  // parallel vector and cray routines. Slow and scalar routines can use any
  // positive integer and thus dihedral constrains can be of any periodicity.
  // Reference angle 0.0 and 180.0 degree correspond to minimum in staggered
  // and eclipsed respectively. Any reference angle is allowed. The value
  // 180 should be prefered over -180 since it is parsed faster and more
  // accuratly. When the periodicity is given as zero, for OTHER THAN THE
  // FIRST dihdral in a multiple dihedral set, then a the amplitude is a
  // constant added to the energy. This is needed to effect the
  // Ryckaert-Bellemans potential for hydrocarbons (see below).

  // First create look-up-table
  map<string, vector<PAR::Improper>::const_iterator> improperLookUpTable;
  for (vector<PAR::Improper>::const_iterator improper = par.impropers.begin();
       improper != par.impropers.end(); improper++)
    improperLookUpTable[String::toUpper(improper->atom1) + "," + String::toUpper(improper->atom2) + "," +
                        String::toUpper(improper->atom3) + "," + String::toUpper(improper->atom4)] = improper;

  // Find the parameters from PAR
  // loop over the improper list in the PSF object
  for (vector<PSF::Improper>::const_iterator improper = psf.impropers.begin();
       improper != psf.impropers.end(); improper++) {
    // store the ID numbers of the atoms in this improper
    int atom1 = improper->atom1 - 1;
    int atom2 = improper->atom2 - 1;
    int atom3 = improper->atom3 - 1;
    int atom4 = improper->atom4 - 1;

    // store the type names of the atoms in this improper
    string improper1 = topo->atomTypes[topo->atoms[atom1].type].name;
    string improper2 = topo->atomTypes[topo->atoms[atom2].type].name;
    string improper3 = topo->atomTypes[topo->atoms[atom3].type].name;
    string improper4 = topo->atomTypes[topo->atoms[atom4].type].name;

    map<string, vector<PAR::Improper>::const_iterator>::const_iterator
    currentimproper =
      improperLookUpTable.find(improper1 + "," + improper2 + "," + improper3 +
        "," + improper4);
    if (currentimproper == improperLookUpTable.end())
      currentimproper = improperLookUpTable.find(
        improper4 + "," + improper3 + "," + improper2 + "," + improper1);

    // Try wildcards if necessary
    // 2) A - X - X - B
    if (currentimproper == improperLookUpTable.end()) {
      currentimproper = improperLookUpTable.find(
        improper1 + ",X,X," + improper4);
      if (currentimproper == improperLookUpTable.end())
        currentimproper = improperLookUpTable.find(
          improper4 + ",X,X," + improper1);
    }
    // 3) X - A - B - C
    if (currentimproper == improperLookUpTable.end()) {
      currentimproper = improperLookUpTable.find(
        "X," + improper2 + "," + improper3 + "," + improper4);
      if (currentimproper == improperLookUpTable.end())
        currentimproper = improperLookUpTable.find(
          improper4 + "," + improper3 + "," + improper2 + ",X");
    }

    // 4) X - A - B - X
    if (currentimproper == improperLookUpTable.end()) {
      currentimproper = improperLookUpTable.find(
        "X," + improper2 + "," + improper3 + ",X");
      if (currentimproper == improperLookUpTable.end())
        currentimproper = improperLookUpTable.find(
          "X," + improper3 + "," + improper2 + ",X");
    }

    // 5) X - X - A - B
    if (currentimproper == improperLookUpTable.end()) {
      currentimproper = improperLookUpTable.find(
        "X,X," + improper3 + "," + improper4);
      if (currentimproper == improperLookUpTable.end())
        currentimproper = improperLookUpTable.find(
          improper4 + "," + improper3 + ",X,X");
    }

    // if we still have not found this improper type in the PAR object, report
    // an error
    if (currentimproper == improperLookUpTable.end())
      THROW("Could not find improper.");

    // if we have found this improper type then copy the
    // improper parameters into the topology
    Torsion torsion;
    torsion.atom1 = atom1;
    torsion.atom2 = atom2;
    torsion.atom3 = atom3;
    torsion.atom4 = atom4;
    torsion.periodicity.push_back(currentimproper->second->periodicity);
    torsion.forceConstant.push_back(currentimproper->second->forceConstant);
    torsion.phaseShift.push_back(dtor(currentimproper->second->phaseShift));
    torsion.multiplicity = 1;
    topo->impropers.push_back(torsion);
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Ryckert-Belleman Dihedral
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  map<string, vector<PAR::RBDihedral>::const_iterator> rb_dihedralLookUpTable;

  for (vector<PAR::RBDihedral>::const_iterator rbdihe = par.rb_dihedrals.begin();
       rbdihe != par.rb_dihedrals.end(); ++rbdihe)
    rb_dihedralLookUpTable[rbdihe->atom1 + "," + rbdihe->atom2 + "," +
                        rbdihe->atom3 + "," + rbdihe->atom4] = rbdihe;

  for (vector<PSF::RBDihedral>::const_iterator rb_dihe = psf.rb_dihedrals.begin();
       rb_dihe != psf.rb_dihedrals.end(); ++rb_dihe) {

    // store the ID numbers of the atoms in this dihedral
    int atom1 = rb_dihe->atom1 - 1;
    int atom2 = rb_dihe->atom2 - 1;
    int atom3 = rb_dihe->atom3 - 1;
    int atom4 = rb_dihe->atom4 - 1;

    // store the type names of the atoms in this dihedral
    string dihedral1 = topo->atomTypes[topo->atoms[atom1].type].name;
    string dihedral2 = topo->atomTypes[topo->atoms[atom2].type].name;
    string dihedral3 = topo->atomTypes[topo->atoms[atom3].type].name;
    string dihedral4 = topo->atomTypes[topo->atoms[atom4].type].name;

    map<string,
        vector<PAR::RBDihedral>::const_iterator>::const_iterator
    currentdihedral =
      rb_dihedralLookUpTable.find(
        dihedral1 + "," + dihedral2 + "," + dihedral3 + "," + dihedral4);

    // if this dihedral type has not been found, try reversing the order of
    // the atom types
    if (currentdihedral == rb_dihedralLookUpTable.end())
      currentdihedral = rb_dihedralLookUpTable.find(
        dihedral4 + "," + dihedral3 + "," + dihedral2 + "," + dihedral1);

    // Try wildcards if necessary
    if (currentdihedral == rb_dihedralLookUpTable.end()) {
      currentdihedral = rb_dihedralLookUpTable.find(
        "X," + dihedral2 + "," + dihedral3 + ",X");
      if (currentdihedral == rb_dihedralLookUpTable.end())
        currentdihedral = rb_dihedralLookUpTable.find(
          "X," + dihedral3 + "," + dihedral2 + ",X");

     //for GROMACS
      if (currentdihedral == rb_dihedralLookUpTable.end())
         currentdihedral = rb_dihedralLookUpTable.find(
         "X," + dihedral2 + "," + dihedral3 + "," + dihedral4);

    }

    // if we still have not found this dihedral type in the PAR object, report
    // an error
    if (currentdihedral == rb_dihedralLookUpTable.end())
      THROWS(
        "Could not find dihedral '" << dihedral1 << "'-'" << dihedral2
                                    << "'-'" << dihedral3 << "'-'" <<
        dihedral4 << "'.");

    // if we have found this dihedral type then copy the
    // dihedral parameters into the topology
    RBTorsion rb_torsion;
    rb_torsion.atom1 = atom1;
    rb_torsion.atom2 = atom2;
    rb_torsion.atom3 = atom3;
    rb_torsion.atom4 = atom4;

    rb_torsion.C0  = currentdihedral->second->C0;
    rb_torsion.C1  = currentdihedral->second->C1;
    rb_torsion.C2  = currentdihedral->second->C2;
    rb_torsion.C3  = currentdihedral->second->C3;
    rb_torsion.C4  = currentdihedral->second->C4;
    rb_torsion.C5  = currentdihedral->second->C5;

    topo->rb_dihedrals.push_back(rb_torsion);

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // LennardJonesParameters
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //
  //Add factors for modifying 1-4 interaction
  //
  if (topo->forceFieldFlag == GROMACS) {
    topo->coulombScalingFactor = par.fudgeQQ;
    topo->LJScalingFactor = par.fudgeLJ;
  }
        
  // get some array sizes
  unsigned int sizeAtomTypes = topo->atomTypes.size();
  unsigned int sizeNonbondeds = par.nonbondeds.size();
  unsigned int sizeNbfixs = par.nbfixs.size();

  topo->lennardJonesParameters.resize(sizeAtomTypes);

  // Nonbonded
  for (unsigned int i = 0; i < sizeAtomTypes; i++) {
    int ti = 0;
    unsigned int bi = sizeNonbondeds;

    for (unsigned int k = 0; k < sizeNonbondeds; ++k) {
      int ok = equalWildcard(par.nonbondeds[k].atom, topo->atomTypes[i].name);
      if (ok > ti) {
        bi = k;
        ti = ok;
      }

      if (ti > 1) break;
    }

    if (ti <= 0)
      THROWS("Could not find matching parameter nonbonded of atom '"
            << topo->atomTypes[i].name << "'.");

    //for implicit solvents we require the van der waals radius from the LJ params
    topo->atomTypes[i].vdwR = par.nonbondeds[bi].sigma;
    //Require all parameters for Amber and OpenMM
    topo->atomTypes[i].sigma = par.nonbondeds[bi].sigma;
    topo->atomTypes[i].sigma14 = par.nonbondeds[bi].sigma14;
    topo->atomTypes[i].epsilon = par.nonbondeds[bi].epsilon;
    topo->atomTypes[i].epsilon14 = par.nonbondeds[bi].epsilon14;

    for (unsigned int j = i; j < sizeAtomTypes; j++) {
      int tj = 0;
      unsigned int bj = sizeNonbondeds;
      for (unsigned int k = 0; k < sizeNonbondeds; ++k) {
        int ok = equalWildcard(par.nonbondeds[k].atom,
          topo->atomTypes[j].name);
        if (ok > tj) {
          bj = k;
          tj = ok;
        }

        if (tj > 1) break;
      }

      if (tj <= 0)
        THROWS("Could not find matching parameter nonbonded of atom '"
               << topo->atomTypes[j].name << "'.");

      LennardJonesParameters paramsij;

      // Charmm28
      Real sigma_i = par.nonbondeds[bi].sigma;
      Real sigma_j = par.nonbondeds[bj].sigma;
      Real sigma14_i = par.nonbondeds[bi].sigma14;
      Real sigma14_j = par.nonbondeds[bj].sigma14;

      Real epsilon_i = par.nonbondeds[bi].epsilon;
      Real epsilon_j = par.nonbondeds[bj].epsilon;
      Real epsilon14_i = par.nonbondeds[bi].epsilon14;
      Real epsilon14_j = par.nonbondeds[bj].epsilon14;

      Real r_ij, e_ij, r14_ij, e14_ij;

      switch(topo->forceFieldFlag) {

        case GROMACS :  
                  r_ij = 0.5 * (sigma_i + sigma_j);
                  e_ij = sqrt(epsilon_i * epsilon_j);
                  r14_ij = 0.5 * sigma14_i + sigma14_j;
                  e14_ij = sqrt(epsilon14_i * epsilon14_j);

                  paramsij.A = power<12>(r_ij) * e_ij * 4.0;
                  paramsij.B = power<6>(r_ij) * e_ij * 4.0;
                  ///#### FudgeLJ=0.5 should be read from the parameter files
                  paramsij.A14 = topo->LJScalingFactor * paramsij.A;//power<12>(r14_ij) * e14_ij * 4.0;
                  paramsij.B14 = topo->LJScalingFactor * paramsij.B;//2 * power<6>(r14_ij) * e14_ij * 4.0;
                  break;

        default /*CHARMM*/:  
                  r_ij = sigma_i + sigma_j;
                  e_ij = sqrt(epsilon_i * epsilon_j);
                  r14_ij = sigma14_i + sigma14_j;
                  e14_ij = sqrt(epsilon14_i * epsilon14_j);

                  paramsij.A = power<12>(r_ij) * e_ij;
                  paramsij.B = 2 * power<6>(r_ij) * e_ij;
                  paramsij.A14 = power<12>(r14_ij) * e14_ij;
                  paramsij.B14 = 2 * power<6>(r14_ij) * e14_ij;
                  break;

      }

      topo->lennardJonesParameters.set(i, j, paramsij);
    }
  }


  //check if Gromacs GB data is available
  if (par.gb_parameters.size() > 0) {

      topo->implicitSolvent = GBSA;
  
      unsigned int typesSize = topo->atomTypes.size();

      for (unsigned int i=0;i<typesSize;i++) {
         AtomType *tempatomtype = &topo->atomTypes[i];
         std::string name = tempatomtype->name;

         map<string,PAR::GB_gromacs>::const_iterator gb_iterator = par.gb_parameters.find(name);
         if (gb_iterator != par.gb_parameters.end()) {
             tempatomtype->vdwR = gb_iterator->second.radius;
         }else {
             tempatomtype->vdwR = 0.00001;
             report << hint << "vdwR for GB not found for atom " << name << endr;
         }

      }
      for (unsigned int i=0;i<typesSize;i++) {
         report << debug(2) <<"BuildTopology : atom_type_name "<<topo->atomTypes[i].name<<" vdwR = "<<topo->atomTypes[i].vdwR<<endr;
      }

  }



  // NbFix
  for (unsigned int k = 0; k < sizeNbfixs; ++k) {
    int ti = 0;
    int tj = 0;
    unsigned int bi = sizeNbfixs;
    unsigned int bj = sizeNbfixs;

    for (unsigned int i = 0; i < sizeAtomTypes; i++) {
      int ok = equalWildcard(par.nbfixs[k].atom1, topo->atomTypes[i].name);
      if (ok > ti) {
        bi = i;
        ti = ok;
      }

      if (ti > 1) break;
    }

    if (ti <= 0) continue;

    for (unsigned int j = 0; j < sizeAtomTypes; j++) {
      int ok = equalWildcard(par.nbfixs[k].atom2, topo->atomTypes[j].name);
      if (ok > tj) {
        bj = j;
        tj = ok;
      }

      if (tj > 1) break;
    }

    if (tj <= 0)
      THROWS("Could not find matching parameter nbfix of atoms '"
             << par.nbfixs[k].atom1 << "' - '" << par.nbfixs[k].atom2 << "'.");

    LennardJonesParameters paramsij;

    paramsij.A = par.nbfixs[k].a;
    paramsij.B = par.nbfixs[k].b;
    paramsij.A14 = par.nbfixs[k].a14;
    paramsij.B14 = par.nbfixs[k].b14;
    topo->lennardJonesParameters.set(bi, bj, paramsij);
  }

  // end loop over NbFix types

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // SCPISM data
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (topo->doSCPISM) {

    // Types  
    unsigned int typesSize = topo->atomTypes.size();

    for (unsigned int i = 0; i < typesSize; i++) {    

      AtomType *tempatomtype = &topo->atomTypes[i];
      std::string name = tempatomtype->name;

      if (mySCPISMTable->myData.find(name) !=
        mySCPISMTable->myData.end()) {

        tempatomtype->mySCPISM_T = new SCPISMAtomTypeParameters();

        // fill in data from inputs
        tempatomtype->mySCPISM_T->alpha =
          mySCPISMTable->myData[name].alpha_i;
        tempatomtype->mySCPISM_T->sqrt_alpha =
          sqrt(tempatomtype->mySCPISM_T->alpha);
        tempatomtype->mySCPISM_T->g_i =
          mySCPISMTable->myData[name].hbond_factor;
        tempatomtype->mySCPISM_T->isHbonded =
          mySCPISMTable->myData[name].isHbonded;
        tempatomtype->mySCPISM_T->A_i = mySCPISMTable->myData[name].A_i;
        tempatomtype->mySCPISM_T->B_i = mySCPISMTable->myData[name].B_i;
        tempatomtype->mySCPISM_T->C_i = mySCPISMTable->myData[name].C_i;

      } else {
        THROWS("No SCPISM data for atom type: " << name);

      }

    }

    // Atoms  
    unsigned int atomsSize = topo->atoms.size();
    for (unsigned int i = 0; i < atomsSize; i++) {  
    
      Atom *tempatom = &topo->atoms[i];
      int type = tempatom->type;
      std::string name = topo->atomTypes[type].name;

      tempatom->mySCPISM_A = new SCPISMAtomParameters();

      // calculated values
      Real R_vdw = topo->atomTypes[type].vdwR + 1.40;  //van der waals strored from LJ params
      Real dR_vdw2 = 1.0 / (4.0 * M_PI * R_vdw * R_vdw);
      Real r_cov = mySCPISMTable->myData[name].r_cov;

      //Variables for original implementation, now just for screened Coulombic
      tempatom->mySCPISM_A->sqrtalphaSCPISM = topo->atomTypes[type].mySCPISM_T->sqrt_alpha;

      //Updated method implementation CRS 01/11/09
      tempatom->mySCPISM_A->zeta = r_cov + (tempatom->scaledCharge > 0 ? 0.85 : 0.35) + 0.5 
                                    - dR_vdw2 * 0.5 * topo->atomTypes[type].mySCPISM_T->A_i;
      tempatom->mySCPISM_A->eta = dR_vdw2 * 0.5 * topo->atomTypes[type].mySCPISM_T->B_i;
      tempatom->mySCPISM_A->energySum = false;
      //
      tempatom->mySCPISM_A->bornRadius = 0.0;

    }

  }

  // end of SCPISM

  //GBSA parameters]
  if (topo->doGBSAOpenMM) {
  unsigned int atomsSize = topo->atoms.size();
  //string hname("H"), cname("C"), 
  for (unsigned int i=0;i < atomsSize; i++) {

    Atom *tempatom = &(topo->atoms[i]);

    tempatom->myGBSA_T = new GBSAAtomParameters();

    int type = tempatom->type;

    string name = topo->atomTypes[type].name;

    tempatom->myGBSA_T->offsetRadius = 0.09;
    //tempatom->myGBSA_T->scalingFactor 
    if ( name[0] == 'H') {
       tempatom->myGBSA_T->scalingFactor = 0.85;
    }else if ( name[0] == 'C') {
       tempatom->myGBSA_T->scalingFactor = 0.72;
    }else if ( name[0] == 'N') {
       tempatom->myGBSA_T->scalingFactor = 0.79;
      
    }else if ( name[0] == 'O') {
       tempatom->myGBSA_T->scalingFactor = 0.85;
    } else if ( name[0] == 'P') {
       tempatom->myGBSA_T->scalingFactor = 0.86;
    } else if (name[0] == 'S') {
       tempatom->myGBSA_T->scalingFactor = 0.96;
    } else {
       tempatom->myGBSA_T->scalingFactor = 0.8;
    }

    //allocate the array to store derivatives of born radius w.r.t. r_{ij}'s
    if (tempatom->myGBSA_T->bornRadiusDerivatives == NULL) {
        tempatom->myGBSA_T->SetSpaceForBornRadiusDerivatives(atomsSize);
    }

    if (tempatom->myGBSA_T->Lvalues == NULL) {
        tempatom->myGBSA_T->SetSpaceLvalues(atomsSize);
    }
    if (tempatom->myGBSA_T->Uvalues == NULL) {
        tempatom->myGBSA_T->SetSpaceUvalues(atomsSize);
    }

    if (tempatom->myGBSA_T->distij == NULL) {
       tempatom->myGBSA_T->SetSpaceDistij(atomsSize);
    }

    // Van der Waal Radii can differ for two atoms of the same type. As part
    // of an incremental fix, we've created a new field in Atom for use with GB.
    // We are not removing the field in AtomType until SCPISM has been fixed to
    // use the new field. In the mean time, copy values from AtomType to Atom.
    tempatom->myGBSA_T->vanDerWaalRadius = topo->atomTypes[type].vdwR;

  }

  }
  
  
  

  // store the molecule information
  buildMoleculeTable(topo);
  buildExclusionTable(topo, topo->exclude);
}

//____~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//____
//____  buildMoleculeTable
//____
//____~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void ProtoMol::buildMoleculeTable(GenericTopology *topo) {
  // *** First we clear all molecules ***
  topo->molecules.clear();

  const unsigned int numAtoms = topo->atoms.size();

  // *** Collecting all possible connections, building the graph ***
  vector<vector<int> > graph(numAtoms, vector<int>());
  set<pair<int, int> > pairs;
  // *** Bonds ***
  for (unsigned int i = 0; i < topo->bonds.size(); i++) {
    int a1 = topo->bonds[i].atom1;
    int a2 = topo->bonds[i].atom2;
    graph[a1].push_back(a2);
    graph[a2].push_back(a1);
    pairs.insert(pair<int, int>(min(a1, a2), max(a1, a2)));
  }

  unsigned int count = pairs.size();

  // *** Angles ***
  for (unsigned int i = 0; i < topo->angles.size(); i++) {
    int a1 = topo->angles[i].atom1;
    int a2 = topo->angles[i].atom2;
    int a3 = topo->angles[i].atom3;
    graph[a1].push_back(a2);
    graph[a1].push_back(a3);
    graph[a2].push_back(a1);
    graph[a2].push_back(a3);
    graph[a3].push_back(a1);
    graph[a3].push_back(a2);
    pairs.insert(pair<int, int>(min(a1, a2), max(a1, a2)));
    pairs.insert(pair<int, int>(min(a3, a2), max(a3, a2)));
  }

  if (count < pairs.size())
    report << hint << "Angles added " << pairs.size() - count <<
    " new bond(s)." << endr;
  count = pairs.size();

  // *** Dihedrals ***
  for (unsigned int i = 0; i < topo->dihedrals.size(); i++) {
    int a1 = topo->dihedrals[i].atom1;
    int a2 = topo->dihedrals[i].atom2;
    int a3 = topo->dihedrals[i].atom3;
    int a4 = topo->dihedrals[i].atom4;
    graph[a1].push_back(a2);
    graph[a1].push_back(a3);
    graph[a1].push_back(a4);
    graph[a2].push_back(a1);
    graph[a2].push_back(a3);
    graph[a2].push_back(a4);
    graph[a3].push_back(a1);
    graph[a3].push_back(a2);
    graph[a3].push_back(a4);
    graph[a4].push_back(a1);
    graph[a4].push_back(a2);
    graph[a4].push_back(a3);
    pairs.insert(pair<int, int>(min(a1, a2), max(a1, a2)));
    pairs.insert(pair<int, int>(min(a3, a2), max(a3, a2)));
    pairs.insert(pair<int, int>(min(a3, a4), max(a3, a4)));
  }

  if (count < pairs.size())
    report << hint << "Dihedrals added " << pairs.size() - count <<
    " new bond(s)." << endr;
  count = pairs.size();

  // *** Impropers ***
  set<pair<int, int> > pairsAddImpropers;
  // Impropers are defined over the bonds 1-2,1-3,1-4 or 4-1,4-3,4-2
  // but MTorsionSystemForce computes distances betweeen 1-2,2-3,3-4
  // we have to take care about these differences ...
  for (unsigned int i = 0; i < topo->impropers.size(); i++) {
    int a1 = topo->impropers[i].atom1;
    int a2 = topo->impropers[i].atom2;
    int a3 = topo->impropers[i].atom3;
    int a4 = topo->impropers[i].atom4;
    graph[a1].push_back(a2);
    graph[a1].push_back(a3);
    graph[a1].push_back(a4);
    graph[a2].push_back(a1);
    graph[a2].push_back(a3);
    graph[a2].push_back(a4);
    graph[a3].push_back(a1);
    graph[a3].push_back(a2);
    graph[a3].push_back(a4);
    graph[a4].push_back(a1);
    graph[a4].push_back(a2);
    graph[a4].push_back(a3);
    pair<int, int> p0(min(a1, a2), max(a1, a2));
    pair<int, int> p1(min(a1, a3), max(a1, a3));
    pair<int, int> p2(min(a1, a4), max(a1, a4));
    pair<int, int> p3(min(a2, a3), max(a2, a3));
    pair<int, int> p4(min(a2, a4), max(a2, a4));
    pair<int, int> p5(min(a3, a4), max(a3, a4));
    int j0 = 0;
    int j1 = 0;
    int j2 = 0;
    int j3 = 0;
    int j4 = 0;
    int j5 = 0;
    if (pairs.find(p0) != pairs.end()) j0++;
    if (pairs.find(p1) != pairs.end()) j1++;
    if (pairs.find(p2) != pairs.end()) j2++;
    if (pairs.find(p3) != pairs.end()) j3++;
    if (pairs.find(p4) != pairs.end()) j4++;
    if (pairs.find(p5) != pairs.end()) j5++;
    if (j0 + j1 + j2 + j3 + j4 + j5 < 3) {
      pairs.insert(p0);
      pairs.insert(p1);
      pairs.insert(p2);
    }
    pairsAddImpropers.insert(p0);
    pairsAddImpropers.insert(p3);
    pairsAddImpropers.insert(p5);
  }

  if (count < pairs.size())
    report << hint << "Impropers added " << pairs.size() - count
           << " new bond(s)." << endr;

  // Now add the improper pairs
  for (set<pair<int, int> >::const_iterator i = pairsAddImpropers.begin();
       i != pairsAddImpropers.end(); i++)
    pairs.insert(*i);

  count = pairs.size();
  // To keep track which atoms already have been added
  // to molecules.
  vector<char> unused(numAtoms, 1);

  // Recursively finding the atoms beloning to a molecule
  for (unsigned int i = 0; i < numAtoms; i++) {
    vector<int> v;
    vector<PairInt> p;
    findNextNeighbor(i, v, p, unused, graph, pairs);

    if (!v.empty()) {
      sort(v.begin(), v.end());
      // add this atom list to the molecules array
      Molecule mol;
      mol.atoms = v;
      for (unsigned int j = 0; j < p.size(); j++)
        if (p[j].first > p[j].second)
          swap(p[j].first, p[j].second);

      sort(p.begin(), p.end());
      mol.pairs = p;
      topo->molecules.push_back(mol);
    }
  }

  // Uncomment to sort descending after size()
  // sort(topo->molecules.begin(),topo->molecules.end(),cmpSize);

  // Look up table for atoms
  const string h("H");
  const string o("O");
  for (unsigned int i = 0; i < topo->molecules.size(); i++) {
    Real mass = 0.0;
    const vector<int> &mol = topo->molecules[i].atoms;
    for (unsigned int j = 0; j < mol.size(); j++) {
      int k = mol[j];
      topo->atoms[k].molecule = i;
      mass += topo->atoms[k].scaledMass;
    }

    topo->molecules[i].mass = mass;
    topo->molecules[i].water =
      (mol.size() == 3 &&
       ((topo->atomTypes[topo->atoms[mol[0]].type].symbolName == h &&
         topo->atomTypes[topo->atoms[mol[1]].type].symbolName == h &&
         topo->atomTypes[topo->atoms[mol[2]].type].symbolName == o) ||
        (topo->atomTypes[topo->atoms[mol[0]].type].symbolName == h &&
         topo->atomTypes[topo->atoms[mol[1]].type].symbolName == o &&
         topo->atomTypes[topo->atoms[mol[2]].type].symbolName == h) ||
        (topo->atomTypes[topo->atoms[mol[0]].type].symbolName == o &&
         topo->atomTypes[topo->atoms[mol[1]].type].symbolName == h &&
         topo->atomTypes[topo->atoms[mol[2]].type].symbolName == h)));
  }

#if defined (DEBUG_PRINT_MOLECULETABLE)
  report << plain << endl
         << "[buildMoleculeTable]: molecule table printout:" << endl;

  for (int i = 0; i < topo->molecules.size(); i++) {
    for (int j = 0; j < topo->molecules[i].size(); j++)
      report << topo->molecules[i][j] << " ";

    report << endl;
  }

  report << endr;
#endif
}

//____~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//____
//____  buildExclusionTable
//____
//____~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void ProtoMol::buildExclusionTable(GenericTopology *topo,
                                   const ExclusionType &exclusionType) {
  if (!exclusionType.valid())
    THROW("[buildExclusionTable()] Exclusion type not defined/valid.");

  topo->exclude = exclusionType;

  //  Resize array.
  topo->exclusions.resize(topo->atoms.size());

  //  If exclusionType is equal to NONE, return.
  if (exclusionType == ExclusionType::NONE) return;

  const int numBonds = topo->bonds.size(),
    numAngles = topo->angles.size(),
    numDihedrals = topo->dihedrals.size(),
    numRBDihedrals = topo->rb_dihedrals.size();

  //  Add excluded bonds.
  for (int i = 0; i < numBonds; i++)
    topo->exclusions.add(topo->bonds[i].atom1, topo->bonds[i].atom2,
      EXCLUSION_FULL);

  if (exclusionType != ExclusionType::ONE2) {
    //  Add excluded angles.
    for (int i = 0; i < numAngles; i++)
      topo->exclusions.add(topo->angles[i].atom1,
        topo->angles[i].atom3, EXCLUSION_FULL);

    if (exclusionType != ExclusionType::ONE3) {

      //  Add excluded dihedrals.
      for (int i = 0; i < numDihedrals; i++) {

        if (exclusionType == ExclusionType::ONE4)
          topo->exclusions.add(topo->dihedrals[i].atom1,
            topo->dihedrals[i].atom4, EXCLUSION_FULL);
        else
          topo->exclusions.add(topo->dihedrals[i].atom1,
            topo->dihedrals[i].atom4, EXCLUSION_MODIFIED);
      }

      //  Add excluded RB dihedrals.
      for (int i = 0; i < numRBDihedrals; i++) {

        if (exclusionType == ExclusionType::ONE4)
          topo->exclusions.add(topo->rb_dihedrals[i].atom1,
            topo->rb_dihedrals[i].atom4, EXCLUSION_FULL);
        else
          topo->exclusions.add(topo->rb_dihedrals[i].atom1,
            topo->rb_dihedrals[i].atom4, EXCLUSION_MODIFIED);
      }
    }

  }

  topo->exclusions.optimize();
}
