#ifndef PROTOMOL_OUTPUT_XTC_TRAJECTORY_H
#define PROTOMOL_OUTPUT_XTC_TRAJECTORY_H

#include "Output.h"

namespace ProtoMol {
  class OutputXTCTrajectory : public Output {
    void *fxtc;
    bool minimalImage;
    int frameOffset;
    string filename;

  public:
    static const std::string keyword;

    OutputXTCTrajectory();
    OutputXTCTrajectory(const std::string &filename, int freq, bool minimal,
                        int frameoffs);
    virtual ~OutputXTCTrajectory() {}

    //  From Output
    Output *doMake(const std::vector<Value> &values) const;
    void doInitialize();
    void doRun(int step);
    void doFinalize(int step);

    //  From Makeabl
    std::string getIdNoAlias() const {return keyword;}
    void getParameters(std::vector<Parameter> &parameter) const;
    bool adjustWithDefaultParameters(std::vector<Value> &values,
                                     const Configuration *config) const;
  };
}
#endif //  PROTOMOL_OUTPUT_XTC_TRAJECTORY_H
