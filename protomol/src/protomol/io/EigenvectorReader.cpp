#include <protomol/io/EigenvectorReader.h>

#include <protomol/base/StringUtilities.h>
#include <protomol/base/Report.h>
#include <protomol/base/MathUtilities.h>
#include <protomol/base/SystemUtilities.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____EigenvectorReader
EigenvectorReader::EigenvectorReader() :
  Reader(ios::in | ios::binary), myEigenvectorInfo(NULL), mySwapEndian(false) {}

EigenvectorReader::EigenvectorReader(const string &filename) :
  Reader(ios::in | ios::binary,
         filename), myEigenvectorInfo(NULL), mySwapEndian(false) {}

EigenvectorReader::EigenvectorReader(const char *filename) :
  Reader(ios::in | ios::binary,
         string(filename)), myEigenvectorInfo(NULL), mySwapEndian(false) {}

EigenvectorReader::~EigenvectorReader() {
  if (myEigenvectorInfo != NULL)
    delete myEigenvectorInfo;
}

bool EigenvectorReader::read() {
  if (myEigenvectorInfo == NULL)
    myEigenvectorInfo = new EigenvectorInfo();
  return read(*myEigenvectorInfo);
}

bool EigenvectorReader::read(EigenvectorInfo &ei) {
  if (!open())
    return false;

  int32 el, ne;
  File::read(reinterpret_cast<char *>(&el), sizeof(int32));
  if (ISLITTLEENDIAN) swapBytes(el);
  ei.myEigenvectorLength = el;
  File::read(reinterpret_cast<char *>(&ne), sizeof(int32));
  if (ISLITTLEENDIAN) swapBytes(ne);
  ei.myNumEigenvectors = ne;
  File::read(reinterpret_cast<char *>(&ei.myMaxEigenvalue), sizeof(double));
  if (ISLITTLEENDIAN) swapBytes(ei.myMaxEigenvalue);
  if(!ei.initializeEigenvectors()) return false;
  for (unsigned int i = 0; i < ei.myEigenvectorLength; i++)
    for (unsigned int j = 0; j < ei.myNumEigenvectors; j++)
      for (unsigned int k = 0; k < 3; k++) {   // X, Y, Z
        file.read(
          reinterpret_cast<char *>(&(ei.myEigenvectors[i * 3 *
                                                       ei.
                                                         myNumEigenvectors
                                                       + j * 3 + k])),
          sizeof(double));
        if (ISLITTLEENDIAN) swapBytes(ei.myEigenvectors[i * 3 *
                                                        ei.myNumEigenvectors
                                                        + j * 3 + k]);
      }

  return true;
}

namespace ProtoMol {
  EigenvectorReader &operator>>(EigenvectorReader &eigenvectorReader,
                                EigenvectorInfo &info) {
    eigenvectorReader.read(info);
    return eigenvectorReader;
  }
}

