/*  -*- c++ -*-  */
#ifndef PROTOMOL_OUTPUT_FINAL_PDB_POS_H
#define PROTOMOL_OUTPUT_FINAL_PDB_POS_H

#include "Output.h"

namespace ProtoMol {
  class Configuration;

  class OutputFinalPDBPos : public Output {
    static const std::string keyword;

  protected:
    std::string filename;
    bool minimalImage;

  public:
    OutputFinalPDBPos();
    OutputFinalPDBPos(const std::string &filename, bool minimal);

    //   From class Output
  private:
    Output *doMake(const std::vector<Value> &values) const;
    void doInitialize() {};
    void doRun(int) {};
    void doFinalize(int);

    //  From class Makeabl
  public:
    std::string getIdNoAlias() const {return keyword;}
    void getParameters(std::vector<Parameter> &parameter) const;
    bool adjustWithDefaultParameters(std::vector<Value> &values,
                                     const Configuration *config) const;
  };
}
#endif //  PROTOMOL_OUTPUT_FINAL_PDB_POS_H
