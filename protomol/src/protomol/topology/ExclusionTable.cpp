#include <protomol/topology/ExclusionTable.h>
#include <protomol/base/Report.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ ExclusionPair
ExclusionPair::ExclusionPair() :
  a1(-1), a2(-1), excl(EXCLUSION_NONE) {}

ExclusionPair::ExclusionPair(int a,
                             int b) :
  a1(a), a2(b), excl(EXCLUSION_NONE) {}

ExclusionPair::ExclusionPair(int a, int b,
                             ExclusionClass c) :
  a1(a), a2(b), excl(c) {}

bool ExclusionPair::operator<(const ExclusionPair &e) const {
  if (a1 < e.a1)
    return true;
  else if (a1 > e.a1)
    return false;
  else
    return a2 < e.a2;
}

//____ ExclusionTable

ExclusionTable::ExclusionTable() :
  lowDeltas(NULL), anyExclsForThisDelta(NULL), myMaxDelta(0), myCurrentSize(0),
  myFastDeltaMax(0) {}

ExclusionTable::~ExclusionTable() {
  if (lowDeltas != NULL)
    delete[] lowDeltas;
  if (anyExclsForThisDelta != NULL)
    delete[] anyExclsForThisDelta;
}

void ExclusionTable::resize(int count) {
  if (count < 0)
    count = 0;

  if (count == myCurrentSize) {
    clear();
    return;
  }

  // release
  if (lowDeltas != NULL)
    delete[] lowDeltas;
  lowDeltas = NULL;
  if (anyExclsForThisDelta != NULL)
    delete[] anyExclsForThisDelta;
  anyExclsForThisDelta = NULL;

  if (count > 0) {
    // allocate
    myFastDeltaMax = Constant::FASTDELTAMAX;
    lowDeltas = new ExclusionClass[count * myFastDeltaMax];
    anyExclsForThisDelta = new char[count];
  }

  myCurrentSize = count;
  clear();
}

void ExclusionTable::optimize() {
  mySet.clear();
  int count = myMaxDelta + 1;
  // minimize the lowDeltas array
  if (count < myFastDeltaMax) {
    ExclusionClass *tmp = new ExclusionClass[count * myCurrentSize];
    for (int i = 0; i < myCurrentSize; i++)
      for (int j = 0; j < count; j++)
        tmp[i * count + j] = lowDeltas[i * myFastDeltaMax + j];

    delete[] lowDeltas;
    lowDeltas = tmp;
    report << hint << "Reduced fast delta max of ExclusionTable from " <<
    myFastDeltaMax - 1 << " to " << myMaxDelta << "." << endr;
    myFastDeltaMax = count;
  } else
    report << hint << "Couldn't reduce fast delta max " << myMaxDelta <<
    " of ExclusionTable (" << myFastDeltaMax - 1 << ")." << endr;
}

void ExclusionTable::clear(void) {
  for (int i = 0; i < myCurrentSize * myFastDeltaMax; i++)
    lowDeltas[i] = EXCLUSION_NONE;

  highDeltas.clear();
  for (int i = 0; i < myCurrentSize; i++)
    anyExclsForThisDelta[i] = 0;

  myMaxDelta = 0;
  myTable.resize(0);
  mySet.clear();
}

void ExclusionTable::add(int atom1, int atom2, ExclusionClass type) {
  if (atom1 > atom2)
    swap(atom1, atom2);
  int delta = atom2 - atom1;
  if (type != EXCLUSION_NONE && atom1 != atom2 &&
      mySet.find(ExclusionPair(atom1, atom2, type)) == mySet.end()) {
    myTable.push_back(ExclusionPair(atom1, atom2, type));
    mySet.insert(ExclusionPair(atom1, atom2, type));
  }
  if (delta > myMaxDelta)
    myMaxDelta = delta;
  if (delta < myFastDeltaMax)
    lowDeltas[atom1 * myFastDeltaMax + delta] = type;
  else {
    highDeltas[PairInt(atom1, delta)] = type;
    anyExclsForThisDelta[delta] = 1;
  }
}

