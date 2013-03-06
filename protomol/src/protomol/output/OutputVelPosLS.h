#ifndef _OUTPUT_VEL_POS_LS_H
#define _OUTPUT_VEL_POS_LS_H

#include <protomol/output/OutputVelPos.h>

namespace ProtoMol {
  using namespace std;

  class OutputVelPosLS : public OutputVelPos
  {
  public:
    static const std::string keyword;

    OutputVelPosLS();
    ~OutputVelPosLS();
    OutputVelPosLS(const string& filename);

    
    string getIdNoAlias() const {return keyword;};

  };
}

#endif
