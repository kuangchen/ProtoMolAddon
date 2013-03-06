#include <protomol/io/DCDTrajectoryWriter.h>

#include <protomol/base/Report.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/base/Exception.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____DCDTrajectoryWriter

DCDTrajectoryWriter::DCDTrajectoryWriter(Real timestep, unsigned int firststep,
                                         bool isLittleEndian) :
  Writer(ios::binary | ios::trunc), myFrameOffset(0),
  myIsLittleEndian(isLittleEndian),
  myFirstStep(firststep), myTimeStep(timestep), 
  firstWrite(true) {}

DCDTrajectoryWriter::DCDTrajectoryWriter(const string &filename, Real timestep,
                                         unsigned int firststep,
                                         bool isLittleEndian) :
  Writer(ios::binary | ios::trunc, filename), myFrameOffset(0),
  myIsLittleEndian(isLittleEndian), myFirstStep(firststep),
  myTimeStep(timestep), firstWrite(true) {}

DCDTrajectoryWriter::DCDTrajectoryWriter(std::ios::openmode mode, int frameoffs,
                                         const string &filename, Real timestep,
                                         unsigned int firststep,
                                         bool isLittleEndian) :
  Writer(mode, filename), myFrameOffset(frameoffs),
  myIsLittleEndian(isLittleEndian), myFirstStep(firststep),
  myTimeStep(timestep), firstWrite(true) {}

DCDTrajectoryWriter::DCDTrajectoryWriter(const char *filename, Real timestep,
                                         unsigned int firststep,
                                         bool isLittleEndian) :
  Writer(ios::binary | ios::trunc, string(filename)), myFrameOffset(0),
  myIsLittleEndian(isLittleEndian), myFirstStep(firststep),
  myTimeStep(timestep), firstWrite(true) {}

bool DCDTrajectoryWriter::openWith(Real timestep, unsigned int firststep,
                                   bool isLittleEndian) {
  setTimestep(timestep);
  setFirststep(firststep);
  setLittleEndian(isLittleEndian);
  return open();
}

bool DCDTrajectoryWriter::openWith(const string &filename, Real timestep,
                                   unsigned int firststep,
                                   bool isLittleEndian) {
  setTimestep(timestep);
  setFirststep(firststep);
  setLittleEndian(isLittleEndian);
  return open(filename);
}

bool DCDTrajectoryWriter::openWith(const char *filename, Real timestep,
                                   unsigned int firststep,
                                   bool isLittleEndian) {
  setTimestep(timestep);
  setFirststep(firststep);
  setLittleEndian(isLittleEndian);
  return open(filename);
}

bool DCDTrajectoryWriter::reopen(unsigned int numAtoms) {
  if (!is_open()) {
    open(filename, ios::binary | ios::in | ios::out);
    if (!is_open()) THROWS("Failed to open '" << filename << "'");
  }

  // Try to read the file size
  file.clear();
  file.seekg(0, ios::end);
  ios::pos_type size = file.tellg();
  if (file.fail()) return false;

  if (size > static_cast<ios::pos_type>(100)) {
    // This is not the first write.  Just update header.

    // 8: Read numSets
    file.seekg(8, ios::beg);
    int32 numSets;
    read((char *)&numSets, 4);
    
    if (myIsLittleEndian != ISLITTLEENDIAN) swapBytes(numSets);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //if frame offset is set, check file length against
    //what we expect. Need to do this each time as file may have
    //numerous additional frames
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //offset for first "end of file"
    streamoff fileoffset = 0;

    //if non zero offset (i.e. checkpoint re-start)
    if( myFrameOffset != 0 ){

      //calculate header size
      const int csize = comment.size() + 9; //9 for "Remarks:"
      int lines = (csize + 80) / 80;        // was 2, now allows for >2 lines
      if((csize % 80) != 0) lines++;
      const int headerSize = 116 + lines * 80;

      //calculate the expected file size
      const unsigned int calculateSize = ( numAtoms * 4 * 3 + 6 * 4 ) * myFrameOffset
                                        + headerSize;

      //report if debug set
      report << debug(2) <<"File size " << size << ", calculated " <<
                calculateSize << ", frame offset " << myFrameOffset << "." << endr;

      //set offset
      fileoffset = calculateSize - size;

      //cannot seek past end, so must be <= 0
      if( fileoffset > 0 )
          THROWS("Corrupt DCD file. Size is " << size << ", should be >= " << calculateSize << ".");

      //correct numsets
      numSets = ++myFrameOffset;
      
    }else{
        //original code here
        ++numSets;
    }

    //back to original code
    
    if (myIsLittleEndian != ISLITTLEENDIAN) swapBytes(numSets);

    //  8: Number of sets of coordinates, NAMD=0 ???
    file.seekp(8, ios::beg);
    file.write((char *)&numSets, 4);  

    // 20: Number of sets of coordinates, NAMD=0 ???
    file.seekp(20, ios::beg);
    file.write((char *)&numSets, 4);   

    file.seekg(fileoffset, ios::end);
    
  } else {
    // First time ...
    file.seekg(0, ios::beg);

    int32 nAtoms = static_cast<int32>(numAtoms);
    int32 numSets = 1;
    int32 numSteps = 1;
    int32 firstStep = static_cast<int32>(myFirstStep);
    float4 timeStep =
      static_cast<float4>(myTimeStep) * Constant::INV_TIMEFACTOR;

    //find comment size
    const int csize = comment.size() + 9; //9 for "Remarks:"

    int32 n0 = 0;
    //int32 n2 = 2;
    int32 n2 = (csize + 80) / 80; // was 2, now allows for >2 lines
    if((csize % 80) != 0) n2++;

    int32 n4 = 4;
    int32 n24 = 24;
    int32 n84 = 84;
    //int32 n164 = 164;
    int32 n164 = n2 * 80 + 4; //was 164, allows for longer comments;

    if (myIsLittleEndian != ISLITTLEENDIAN) {
      swapBytes(nAtoms);
      swapBytes(numSets);
      swapBytes(numSteps);
      swapBytes(firstStep);
      swapBytes(timeStep);

      swapBytes(n0);
      swapBytes(n2);
      swapBytes(n4);
      swapBytes(n24);
      swapBytes(n84);
      swapBytes(n164);
    }

    // Write header
    file.write((char *)&n84, 4); //  0
    file.write(string("CORD").c_str(), 4); //  4
    //  8: Number of sets of coordinates, NAMD=0 ???
    file.write((char *)&numSets, 4);
    // 12: Starting timestep of DCD file, should never be zero
    file.write((char *)&firstStep, 4);
    // 16: Timesteps between DCD saves 
    file.write((char *)&numSteps, 4);
    // 20: NAMD writes += numSteps 
    file.write((char *)&numSets, 4);
    // 24
    for (unsigned i = 0; i < 5; i++)
      file.write((char *)&n0, 4);
    // 44 : length of a timestep
    file.write((char *)&timeStep, 4);
    // 48 : unit cell, none=0, used=1
    for (unsigned i = 0; i < 9; i++)
      file.write((char *)&n0, 4);
    // 84: Pretend to be Charmm 24
    file.write((char *)&n24, 4);   
    file.write((char *)&n84, 4);

    // 92: Write DCD title record
    file.write((char *)&n164, 4);
    file.write((char *)&n2, 4);
    string remarks = string("Remarks: File '") + filename + "'. ProtoMol (" +
      __DATE__ + " at " + __TIME__ + ")";
    file.write(getRightFill(remarks, 80).c_str(), 80);

    //file.write(getRightFill(string("Remarks: " + comment), 80).c_str(), 80);
    const int numChars = (n2 - 1) * 80; //alow longer comments
    file.write(getRightFill(string("Remarks: " + comment), numChars).c_str(), numChars);

    file.write((char *)&n164, 4);

    // 264: Write DCD num-atoms record
    file.write((char *)&n4, 4);
    file.write((char *)&nAtoms, 4);
    file.write((char *)&n4, 4);
  }

  return !file.fail();
}

bool DCDTrajectoryWriter::write(const Vector3DBlock &coords) {

  //don't write first frame if checkpoint re-start
  if(firstWrite && myFrameOffset != 0 ){
    //one shot
    firstWrite = false;
    return true;
  }

  //original code
  const unsigned int count = coords.size();
  if (!reopen(count)) return false;

  myX.resize(count);
  myY.resize(count);
  myZ.resize(count);

  for (unsigned int i = 0; i < count; ++i) {
    myX[i] = static_cast<float>(coords[i].c[0]);
    myY[i] = static_cast<float>(coords[i].c[1]);
    myZ[i] = static_cast<float>(coords[i].c[2]);
    if (myIsLittleEndian != ISLITTLEENDIAN) {
      swapBytes(myX[i]);
      swapBytes(myY[i]);
      swapBytes(myZ[i]);
    }
  }

  int32 nAtoms = static_cast<int32>(count * 4);
  if (myIsLittleEndian != ISLITTLEENDIAN) swapBytes(nAtoms);

  file.write((char *)&nAtoms, sizeof(int32));
  file.write((char *)&(myX[0]), count * sizeof(float4));
  file.write((char *)&nAtoms, sizeof(int32));

  file.write((char *)&nAtoms, sizeof(int32));
  file.write((char *)&(myY[0]), count * sizeof(float4));
  file.write((char *)&nAtoms, sizeof(int32));

  file.write((char *)&nAtoms, sizeof(int32));
  file.write((char *)&(myZ[0]), count * sizeof(float4));
  file.write((char *)&nAtoms, sizeof(int32));

  close();
  return !file.fail();
}

void DCDTrajectoryWriter::setLittleEndian(bool littleEndian) {
  myIsLittleEndian = littleEndian;
}

void DCDTrajectoryWriter::setTimestep(Real timestep) {
  myTimeStep = timestep;
}

void DCDTrajectoryWriter::setFirststep(unsigned int firststep) {
  myFirstStep = firststep;
}


namespace ProtoMol {
  DCDTrajectoryWriter &operator<<(DCDTrajectoryWriter &dcdWriter,
                                  const Vector3DBlock &coords) {
    dcdWriter.write(coords);
    return dcdWriter;
  }

  DCDTrajectoryWriter &operator<<(DCDTrajectoryWriter &dcdWriter,
                                  const XYZ &xyz) {
    dcdWriter.write(xyz.coords);
    return dcdWriter;
  }
}
