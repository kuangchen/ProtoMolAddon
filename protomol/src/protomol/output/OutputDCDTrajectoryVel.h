/*  -*- c++ -*-  */
#ifndef PROTOMOL_OUTPUT_DCD_TRAJECTORY_VEL_H
#define PROTOMOL_OUTPUT_DCD_TRAJECTORY_VEL_H

#include <protomol/output/Output.h>

namespace ProtoMol {
  class DCDTrajectoryWriter;

  class OutputDCDTrajectoryVel : public Output {
  public:
    static const std::string keyword;

  private:
    DCDTrajectoryWriter *dCD;
    bool minimalImage;

  public:
    OutputDCDTrajectoryVel();
    OutputDCDTrajectoryVel(const std::string &filename, int freq, bool minimal);
    virtual ~OutputDCDTrajectoryVel();

    // From class Output
  private:
    Output *doMake(const std::vector<Value> &values) const;
    void doInitialize();
    void doRun(int step);
    void doFinalize(int step);

    // From class Makeable
  public:
    std::string getIdNoAlias() const {return keyword;}
    void getParameters(std::vector<Parameter> &parameter) const;
    bool adjustWithDefaultParameters(std::vector<Value> &values,
                                     const Configuration *config) const;
  };
}
#endif // PROTOMOL_OUTPUT_DCD_TRAJECTORY_VEL_H
