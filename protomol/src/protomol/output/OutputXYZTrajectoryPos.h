/*  -*- c++ -*-  */
#ifndef PROTOMOL_OUTPUT_XYZ_TRAJECTORY_POS_H
#define PROTOMOL_OUTPUT_XYZ_TRAJECTORY_POS_H

#include <protomol/output/Output.h>

namespace ProtoMol {
  class XYZTrajectoryWriter;

  class OutputXYZTrajectoryPos : public Output {
  public:
    static const std::string keyword;

  private:
    XYZTrajectoryWriter *xYZ;
    bool minimalImage;

  public:
    OutputXYZTrajectoryPos();
    OutputXYZTrajectoryPos(const std::string &filename, int freq, bool minimal);
    virtual ~OutputXYZTrajectoryPos();

    //   From class Output
  private:
    Output *doMake(const std::vector<Value> &values) const;
    void doInitialize();
    void doRun(int step);
    void doFinalize(int step);

    //  From class Makeabl
  public:
    std::string getIdNoAlias() const {return keyword;}
    void getParameters(std::vector<Parameter> &parameter) const;
    bool adjustWithDefaultParameters(std::vector<Value> &values,
                                     const Configuration *config) const;
  };
}
#endif //  PROTOMOL_OUTPUT_XYZ_TRAJECTORY_POS_H
