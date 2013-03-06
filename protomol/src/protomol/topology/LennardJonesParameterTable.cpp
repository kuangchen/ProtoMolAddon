#include <protomol/topology/LennardJonesParameterTable.h>
#include <protomol/base/Report.h>
using namespace ProtoMol::Report;

using namespace ProtoMol;
//____ LennardJonesParameterTable
LennardJonesParameterTable::LennardJonesParameterTable() :
  myCurrentSize(0), myData(NULL) {}

LennardJonesParameterTable::~LennardJonesParameterTable() {
  if (myData != NULL)
    delete[] myData;
}

void LennardJonesParameterTable::resize(int count) {
  if (count < 0)
    count = 0;
  if (count == myCurrentSize)
    return;

  if (myData != NULL)
    delete[] myData;
  myData = NULL;

  if (count > 0)
    myData = new LennardJonesParameters[count * count];

  myCurrentSize = count;
}

void LennardJonesParameterTable::set(int type1, int type2,
                                     const LennardJonesParameters &params) {
  if (type1 < 0 || type2 < 0 || type1 >= myCurrentSize || type2 >=
      myCurrentSize) {
    report << recoverable <<
    "[LennardJonesParameters::set] index of out range!" << endr;
    return;
  }
  myData[type1 * myCurrentSize + type2] = params;
  myData[type2 * myCurrentSize + type1] = params;
}

