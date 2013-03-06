/*  -*- c++ -*-  */
#ifndef PROTOMOL_OUTPUT_FINAL_XYZ_VEL_H
#define PROTOMOL_OUTPUT_FINAL_XYZ_VEL_H

#include <protomol/output/Output.h>

namespace ProtoMol {
  class Configuration;

  class OutputFinalXYZVel : public Output {
  public:
    static const std::string keyword;

  private:
    std::string filename;

  public:
    OutputFinalXYZVel();
    OutputFinalXYZVel(const std::string &filename);

    //   From class Output
  private:
    Output *doMake(const std::vector<Value> &values) const;
    void doInitialize() {};
    void doRun(int) {};
    void doFinalize(int step);

    //  From class Makeabl
  public:
    std::string getIdNoAlias() const {return keyword;}
    void getParameters(std::vector<Parameter> &parameter) const;
  };
}
#endif //  PROTOMOL_OUTPUT_FINAL_XYZ_VEL_H
