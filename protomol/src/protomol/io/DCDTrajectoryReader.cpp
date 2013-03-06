#include <protomol/io/DCDTrajectoryReader.h>

#include <protomol/base/Report.h>
#include <protomol/base/SystemUtilities.h>
#include <protomol/base/Exception.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

DCDTrajectoryReader::DCDTrajectoryReader() : Reader(ios::binary), xyz(0), first(true) {}


DCDTrajectoryReader::DCDTrajectoryReader(const string &filename) :
  Reader(ios::binary, filename), xyz(0), first(true) {}


DCDTrajectoryReader::~DCDTrajectoryReader() {
  if (xyz) delete xyz;
}


bool DCDTrajectoryReader::tryFormat() {
  if (!open()) return false;

  // Get file size
  file.seekg(0, ios::end);
  ios::pos_type size = file.tellg();
  file.seekg(0, ios::beg);


  // Get endian
  int32 n = 0;
  File::read((char *)&n, sizeof(int32));
  int32 m = n;
  swapBytes(m);

  //get header
  swap = false;
  if (m == 84)
    swap = true;

  //back to start for read
  file.seekg(0, ios::beg);

  fortranRead((char *)&header, sizeof(header), "header");
  if (swap) {
    swapBytes(header.frames);
    swapBytes(header.firststep);
    swapBytes(header.freeIndexes);
  }

  // Read DCD magic
  char coord[5];
  //File::read(coord, 4);
  for(int i=0; i<4; i++)
      coord[i] = header.cord[i];
  coord[4] = '\0';
  close();

  // Check it
  if (size >= 104 && (n == 84 || m == 84) && string(coord) == "CORD")
    return !file.fail();


  return false;
}


bool DCDTrajectoryReader::read() {
  if (!xyz) xyz = new Vector3DBlock();
  return read(*xyz);
}




bool DCDTrajectoryReader::read(Vector3DBlock &xyz) {
  try {
    doRead(xyz);

    return true;
  } catch (const Exception &e) {}

  return false;
}


void DCDTrajectoryReader::doRead(Vector3DBlock &xyz) {
  try {
    if (first) {
      if (!is_open()) if (!open()) THROW("Open failed");
      // Check endian
      int32 n = 0;
      File::read((char *)&n, sizeof(int32));
      int32 m = n;
      swapBytes(m);
      if (n != 84 && m != 84) THROW("Invalid DCD header");
      if (m == 84) {
	swap = true;
	report
	  << hint << "[DCDTrajectoryReader::read] Reading "
	  << (ISLITTLEENDIAN ? "big" : "little") << "endian input on "
	  << (ISLITTLEENDIAN ? "little" : "big") << "endian machine." << endr;

      } else swap = false;

      // Get file size
      file.seekg(0, ios::end);
      ios::pos_type filesize = file.tellg();
      file.seekg(0, ios::beg);
      if (filesize < 104) THROWS("Invalid file size = " << filesize);


      fortranRead((char *)&header, sizeof(header), "header");
      if (swap) {
	swapBytes(header.frames);
	swapBytes(header.firststep);
	swapBytes(header.freeIndexes);
      }


      // Check header
      if (string(header.cord, 4) != "CORD") THROW("Header invalid");


      // Read comment
      unsigned int size = 0;
      char *commentData = fortranReadX(0, size, "comment");
      if (size) {
	comment.resize(size + 1);
	memcpy(&comment[0], commentData, size);
	comment[size] = 0;
      }


      // Read number of atoms
      natoms = 0;
      fortranRead((char *)&natoms, sizeof(int32), "# atoms");
      if (swap) swapBytes(natoms);
    }

    // Skip free indexes
    if (header.freeIndexes > 0)
      file.seekg(4 * (natoms - header.freeIndexes + 2), ios::cur);

    if (file.fail())
      THROWS("Error skipping " << header.freeIndexes << "free indexes");


    // Read next frame
    std::vector<float4> x(natoms);
    std::vector<float4> y(natoms);
    std::vector<float4> z(natoms);
    //float4 x, y, z;

    xyz.resize(natoms);
    //xyz.resize(header.frames);

    //    for (int i = 0; i < header.frames; i++) {
    fortranRead((char *)&x[0], natoms * 4, "X dimension");
    fortranRead((char *)&y[0], natoms * 4, "Y dimension");
    fortranRead((char *)&z[0], natoms * 4, "Z dimension");
    for (int j = 0; j < natoms; j++) {
      if (swap) {
	swapBytes(x[j]);
	swapBytes(y[j]);
	swapBytes(z[j]);
      }
      xyz.c[3*j] = x[j];
      xyz.c[3*j+1] = y[j];
      xyz.c[3*j+2] = z[j];
    }
    //}
  } catch (const Exception &e) {
    THROWSC("DCD read error at " << file.tellg() << " in "
            << getFilename(), e);
  }
  first = false;
}

Vector3DBlock *DCDTrajectoryReader::orphanXYZ() {
  Vector3DBlock *tmp = xyz;
  xyz = 0;
  return tmp;
}

void DCDTrajectoryReader::fortranRead(char *data, unsigned int size,
                                      const string &err) {
  unsigned int x = size;
  fortranReadX(data, x, err);
}

char *DCDTrajectoryReader::fortranReadX(char *data, unsigned int &size,
                                        const std::string &err) {
  int32 head;

  // Read record head
  File::read((char *)&head, sizeof(int32));
  if (swap) swapBytes(head);

  if (head <= 0) THROWS("Invalid record header: " << err);

  if (size && size != (unsigned int)head)
    THROWS("Record size " << head << " not what expected " << size << ": "
           << err);

  // Read data
  if (!size) {
    size = (unsigned int)head;
    data = new char[size];
    if (!data) THROWS("Failed to allocate " << size << " bytes");
  }
  File::read(data, size);

  // Read record tail
  int32 tail;
  File::read((char *)&tail, sizeof(int32));
  if (swap) swapBytes(tail);

  if (head != tail)
    THROWS("Record tail " << tail << " does not match head " << head
           << ": " << err);

  if (file.fail()) THROW("IO failure");


  return data;
}

int DCDTrajectoryReader::readFirstStep(){
  return header.firststep;
}

namespace ProtoMol {
  DCDTrajectoryReader &operator>>(DCDTrajectoryReader &reader,
                                  Vector3DBlock &xyz) {
    try {
      reader.doRead(xyz);
      
    } catch (const Exception &e) {
      THROWSC("Failed to read DCD file '" << reader.getFilename() << "'", e);
    }
    
    return reader;
  }
}

