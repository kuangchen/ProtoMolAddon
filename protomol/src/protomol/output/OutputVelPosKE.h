#ifndef _OUTPUT_VEL_POS_KE_H
#define _OUTPUT_VEL_POS_KE_H

#include <protomol/output/OutputVelPos.h>

namespace ProtoMol {
  using namespace std;

  class OutputVelPosKE : public OutputVelPos
  {
  public:
    static const std::string keyword;

    OutputVelPosKE();
    ~OutputVelPosKE();
    OutputVelPosKE(const string& filename);

  public:
 
    string getIdNoAlias() const {return keyword;};
  };
}

#endif
