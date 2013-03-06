#ifndef PROTOMOL_DATA_COMPARATOR_H
#define PROTOMOL_DATA_COMPARATOR_H

#include <string>
#include <vector>
#include <protomol/type/Real.h>
#include <protomol/type/Vector3D.h>

namespace ProtoMol {
  //class Vector3D;
  class Vector3DBlock;
  class XYZ;

  class DataComparator {
  public:
    static Real compare(const Real &data1, const Real &data2);
    //static Real compare(const Vector3D &data1, const Vector3D &data2);
    static Real compare(const Vector3DBlock &data1,
                        const Vector3DBlock &data2,
                        Real tolerance, unsigned int &count);
    static Real compare(const std::vector<XYZ> &data1,
                        const std::vector<XYZ> &data2, Real tolerance,
                        unsigned int &count, unsigned int &divergeFrame);
    static Real compare(const std::string &file1, const std::string &file2,
                        Real tolerance, unsigned int &count,
                        unsigned int &divergeFrame);
    static void read(const std::string &filename, std::vector<XYZ> &data);
  };
}

#endif // PROTOMOL_DATA_COMPARATOR_H

