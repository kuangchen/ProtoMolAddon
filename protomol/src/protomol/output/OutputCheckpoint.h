/*  -*- c++ -*-  */
#ifndef PROTOMOL_OUTPUT_CHECKPIINT_H
#define PROTOMOL_OUTPUT_CHECKPIINT_H

#include <protomol/base/StringUtilities.h>
#include <protomol/output/Output.h>
#include <protomol/base/Timer.h>

namespace ProtoMol {
  class Configuration;

  class OutputCheckpoint : public Output {
  public:
    static const std::string keyword;

  protected:
    int current;
    std::string name;
    std::string posBase, velBase;

  public:
    OutputCheckpoint() : current(0) {}
    OutputCheckpoint(const std::string &name, int freq, int start,
                      const std::string &posbase, const std::string &velbase);

  private:
    void WritePositions(int step);
    void WriteVelocities(int step);
    void WriteConfig(int step);

  public:
    void doIt(int step);

  private:
    //  From class Output
    Output *doMake(const std::vector<Value> &values) const;
    void doInitialize();
    void doRun(int step);
    void doFinalize(int) {}
    bool isIdDefined(const Configuration *config) const;
    bool addDoKeyword() const {return false;}

  public:
    //  From class Makeable
    std::string getIdNoAlias() const {return keyword;}
    void getParameters(std::vector<Parameter> &) const;
    bool adjustWithDefaultParameters(std::vector<Value> &values,
                                     const Configuration *config) const;
  };
}

#endif //  PROTOMOL_OUTPUT_CHECKPIINT_H
