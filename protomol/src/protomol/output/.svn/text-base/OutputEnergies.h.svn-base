/*  -*- c++ -*-  */
#ifndef PROTOMOL_OUTPUT_ENERGIES_H
#define PROTOMOL_OUTPUT_ENERGIES_H

#include "Output.h"

#ifdef HAVE_LIBFAH
#include <cbang/os/File.h>
#else
#include <fstream>
#endif

namespace ProtoMol {
  class OutputEnergies : public Output {
  public:
    static const std::string keyword;

  protected:
#ifdef HAVE_LIBFAH
    cb::File file;
#else
    std::ofstream file;
#endif

    std::string filename;
    bool doMolecularTemperature;

  public:
    OutputEnergies() : doMolecularTemperature(false) {}
    OutputEnergies(const std::string &filename, int freq, bool doMolTemp);

  private:
    //   From class OutputFile
    void doRunCached(int step);

    //   From class Output
    void doInitialize();
    void doRun(int step);
    void doFinalize(int step);

  public:
    //  From class Makeable
    Output *doMake(const std::vector<Value> &values) const;
    std::string getIdNoAlias() const {return keyword;}
    void getParameters(std::vector<Parameter> &parameter) const;
  };
}
#endif //  PROTOMOL_OUTPUT_ENERGIES_H
