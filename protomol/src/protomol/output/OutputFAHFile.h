/*  -*- c++ -*-  */
#ifndef PROTOMOL_OUTPUT_FAH_FILE_H
#define PROTOMOL_OUTPUT_FAH_FILE_H

#include <protomol/output/Output.h>
#include <protomol/io/File.h>
#include <string>

namespace ProtoMol {
  class Configuration;

  class OutputFAHFile : public Output, public File {
  public:
    static const std::string keyword;

  private:
    std::string filename;

  public:
    OutputFAHFile();
    OutputFAHFile(const std::string &filename, int freq);

    //   From class Output
  private:
    Output *doMake(const std::vector<Value> &values) const;
    void doInitialize();
    void doRun(int step);
    void doFinalize(int);
    bool isIdDefined(const Configuration *config) const;
    bool addDoKeyword() const {return false;}

    //  From class Makeabl
  public:
    std::string getIdNoAlias() const {return keyword;}
    void getParameters(std::vector<Parameter> &) const;
  };
}

#endif //  PROTOMOL_OUTPUT_FAH_FILE_H
