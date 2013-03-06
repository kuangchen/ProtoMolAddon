/*  -*- c++ -*-  */
#ifndef PROTOMOL_OUTPUT_XYZ_TRAJECTORY_FORCE_H
#define PROTOMOL_OUTPUT_XYZ_TRAJECTORY_FORCE_H

#include <protomol/output/Output.h>

namespace ProtoMol {
  class XYZTrajectoryWriter;

  class OutputXYZTrajectoryForce : public Output {
  public:
    static const std::string keyword;

  private:
    XYZTrajectoryWriter *xYZ;

  public:
    OutputXYZTrajectoryForce();
    OutputXYZTrajectoryForce(const std::string &filename, int freq);
    virtual ~OutputXYZTrajectoryForce();

  private:
    //   From class Output
    Output *doMake(const std::vector<Value> &values) const;
    void doInitialize();
    void doRun(int step);
    void doFinalize(int step);

    //  From class Makeabl
  public:
    std::string getIdNoAlias() const {return keyword;}
    void getParameters(std::vector<Parameter> &parameter) const;
  };
}
#endif //  PROTOMOL_OUTPUT_XYZ_TRAJECTORY_FORCE_H
