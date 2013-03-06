#include "DataComparator.h"

#include <protomol/base/SystemUtilities.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/base/Exception.h>
#include <protomol/type/String.h>

#include <protomol/io/DCDTrajectoryReader.h>
#include <protomol/io/XYZTrajectoryReader.h>
#include <protomol/io/XYZReader.h>
#include <protomol/io/PDBReader.h>
#include <protomol/type/XYZ.h>

#ifdef DATA_COMPARATOR_STANDALONE
#include <iostream>
#include <iomanip>
#endif

using namespace std;
using namespace ProtoMol;

Real DataComparator::compare(const Real &data1, const Real &data2) {
  return RealAbs(data1 - data2);
}

// Real DataComparator::compare(const Vector3D &data1, const Vector3D &data2) {
//   Real max = 0;
//   Real d;

//   d = compare(data1.c[0], data2.c[0]);
//   if (d > max) max = d;
//   d = compare(data1.c[1], data2.c[1]);
//   if (d > max) max = d;
//   d = compare(data1.c[2], data2.c[2]);
//   if (d > max) max = d;

//   return max;
// }

Real DataComparator::compare(const Vector3DBlock &data1,
                             const Vector3DBlock &data2,
                             Real tolerance, unsigned int &count) {
  Real max = 0;
  Real d;


  if (data1.size() != data2.size())
    THROW("Vector3DBlock sizes don't match");

  for (unsigned int i = 0; i < 3*data1.size(); i++) {
    d = compare(data1.c[i], data2.c[i]);

    if (d > max) max = d;

    if (d > tolerance) count++;
  }

  return max;
}

Real DataComparator::compare(const vector<XYZ> &data1,
                             const vector<XYZ> &data2,
                             Real tolerance, unsigned int &count,
                             unsigned int &divergeFrame) {
  Real max = 0;
  Real d;

  divergeFrame = 0;

  if (data1.size() != data2.size())
    THROWS("Frame sizes not equal " << data1.size() << " != " << data2.size());

  for (unsigned int i = 0; i < data1.size(); i++) {
    d = compare(data1[i].coords, data2[i].coords, tolerance, count);
    if (d > max) max = d;

    if (d > tolerance && !divergeFrame) divergeFrame = i + 1;
  }

  return max;
}

Real DataComparator::compare(const string &file1, const string &file2,
                             Real tolerance, unsigned int &count,
                             unsigned int &divergeFrame) {
  vector<XYZ> data1;
  vector<XYZ> data2;

  // Special case for energies file; no XYZ values
  if (file1.substr(file1.size()-8, 8) == "energies") {
    ifstream infile(file1.c_str(), ios::in);
    ifstream origfile(file2.c_str(), ios::in);
    Real d1, d2, c, max=0;
    while (!infile.eof()) {
      infile >> d1;
      origfile >> d2;
      c = compare(d1, d2);
      if (c > max) max = c;
    }
    return max;
  }
  else {
    read(file1, data1);
    read(file2, data2);
  }
  return compare(data1, data2, tolerance, count, divergeFrame);
}

void DataComparator::read(const string &filename, vector<XYZ> &data) {
  if (!isAccessible(filename)) THROWS("Cannot access '" << filename << "'");

  XYZTrajectoryReader reader(filename);
  data.clear();

  try {
    data.push_back(XYZ());
    reader >> data[0].coords; // Will throw exception if it breaks
    Vector3DBlock coords;
    while (reader.read(coords)) {
      data.push_back(XYZ());
      data[data.size()-1].coords = coords;
    }

  } catch (const Exception &e1) {
    reader.close();
    XYZReader reader(filename);

    try {
      data.clear();
      data.push_back(XYZ());
      reader >> data[data.size()-1];
    } catch (const Exception &e2) {
      reader.close();
      data.pop_back();
      DCDTrajectoryReader reader(filename);

      try {
	data.push_back(XYZ());
	reader >> data[0].coords; // Will throw exception if it breaks
	Vector3DBlock coords;
	while (reader.read(coords)) {
	  data.push_back(XYZ());
	  data[data.size()-1].coords = coords;
	}
	//        reader >> data.coords;
      } catch (const Exception &e3) {

          reader.close();
	  data.clear();
          PDBReader reader(filename);

          try {
	    data.push_back(XYZ());
	    reader >> data[data.size()-1];
          } catch(const Exception &e4) {
        THROWS("Unsupported file format '" << filename << "' due to errors:"
               << endl << e1 << endl << e2 << endl << e3 << endl);
          }
      }
    }
  }

#ifdef DATA_COMPARATOR_STANDALONE
  cout << "Frames: " << setw(5) << data.size();
  if (data.size()) cout << " Atoms: " << setw(5) << data[0].size();
  cout << endl;
#endif // DATA_COMPARATOR_STANDALONE
}

#ifdef DATA_COMPARATOR_STANDALONE
#ifndef _WIN32
#include <protomol/debug/Debugger.h>
#endif

int main(int argc, char *argv[]) {
#ifdef HAVE_STACK_TRACE
  Exception::enableStackTraces = true;
#endif

#ifndef _WIN32
  //Debugger::initStackTrace(argv[0]);
#endif

  bool result = true;

  try {

    if (argc != 3 && argc != 4) {
      cerr << "Syntax: " << argv[0] << " <file1> <file2> [tolerance]" << endl;
      return -1;
    }

    Real tolerance = 0;
    unsigned int count = 0;
    unsigned int divergeFrame = 0;

    if (argc == 4) tolerance = String::parseDouble(argv[3]);

    Real d = DataComparator::compare(argv[1], argv[2], tolerance,
                                     count, divergeFrame);

    if (argc == 4) {
      if (tolerance < d) {
        cout << "Files do not match" << endl
             << "  Tolerance         = " << tolerance << endl
             << "  Max diff          = " << d << endl
             << "  Divergent vectors = " << count << endl
             << "  Divergent frame   = " << divergeFrame << endl;

        result = false;

      } else cout << "Files match" << endl;

    } else {
      cout << "Maximum difference: " << d << endl;
      if (d) result = false;
    }

  } catch (const Exception &e) {
    cerr << e << endl;
    cerr << setw(10) << getFileSize(argv[1]) << " bytes " << argv[1] << endl;
    cerr << setw(10) << getFileSize(argv[2]) << " bytes " << argv[2] << endl;
  }

  return result ? 0 : 1;
}

#endif // DATA_COMPARATOR_STANDALONE
