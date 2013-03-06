#include <protomol/io/Reader.h>

using namespace std;
using namespace ProtoMol;

//____ Reader
Reader::Reader() : File(ios::in) {}
Reader::Reader(const string &filename) : File(ios::in, filename) {}
Reader::Reader(ios::openmode mode) : File(ios::in | mode) {}
Reader::Reader(ios::openmode mode, const string &filename) :
  File(ios::in | mode, filename) {}
