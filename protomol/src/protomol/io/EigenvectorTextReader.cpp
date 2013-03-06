#include <protomol/io/EigenvectorTextReader.h>

#include <protomol/base/StringUtilities.h>
#include <protomol/base/Report.h>
#include <protomol/base/MathUtilities.h>
#include <protomol/base/SystemUtilities.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____EigenvectorTextReader
EigenvectorTextReader::EigenvectorTextReader() :
  Reader(), myEigenvectorInfo(0) {}

EigenvectorTextReader::EigenvectorTextReader(const string &filename) :
  Reader(filename), myEigenvectorInfo(0) {}

EigenvectorTextReader::EigenvectorTextReader(const char *filename) :
  Reader(string(filename)), myEigenvectorInfo(0) {}

EigenvectorTextReader::~EigenvectorTextReader() {
  if (myEigenvectorInfo != 0)
    delete myEigenvectorInfo;
}

bool EigenvectorTextReader::read() {
  if (myEigenvectorInfo == 0)
    myEigenvectorInfo = new EigenvectorInfo();
  return read(*myEigenvectorInfo);
}

bool EigenvectorTextReader::read(EigenvectorInfo &ei) {
  if (!open())
    return false;

  int num, num1;
  double ev;
  file >> num;
  report << plain << num << endr;
  file >> num1;
  report << plain << num1 << endr;

  file >> ev;
  report << plain << ev << endr;

  string str;
  for (unsigned int i = 0; i < 4; i++) {
    file >> str;
    report << plain << str << endr;
  }

  ei.myEigenvectorLength = num;
  ei.myNumEigenvectors = num1;
  ei.myMaxEigenvalue = ev;
  if(!ei.initializeEigenvectors()) return false;

  for (int i = 0; i < num; i++) {
    int x;
    double y;
    for (int j = 0; j < num1; j++) {
      file >> x;

      for (int k = 0; k < 3; k++) {
        file >> y;

        ei.myEigenvectors[i * 3 * num1 + j * 3 + k] = y;
      }
    }
  }

  return true;
}


namespace ProtoMol {
  EigenvectorTextReader &
  operator>>(EigenvectorTextReader &eigenvectorReader, EigenvectorInfo &info) {
    eigenvectorReader.read(info);
    return eigenvectorReader;
  }
}
