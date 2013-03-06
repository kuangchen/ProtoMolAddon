/*  -*- c++ -*-  */
#ifndef PROTOMOL_OUTPUT_FINAL_XYZ_POS_H
#define PROTOMOL_OUTPUT_FINAL_XYZ_POS_H

#include <protomol/output/Output.h>

namespace ProtoMol {
  class Configuration;

  class OutputFinalXYZPos : public Output {
  public:
    static const std::string keyword;

  private:
    std::string filename;
    bool minimalImage;

  public:
    OutputFinalXYZPos();
    OutputFinalXYZPos(const std::string &filename, bool minimal);

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
#endif //  PROTOMOL_OUTPUT_FINAL_XYZ_POS_H
