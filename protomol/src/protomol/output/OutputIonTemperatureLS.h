#ifndef PROTOMOL_OUTPUT_ION_TEMPERATURE_LS_H
#define PROTOMOL_OUTPUT_ION_TEMPERATURE_LS_H

#include <protomol/output/Output.h>

namespace ProtoMol {
  using namespace std;
  class OutputIonTemperatureLS : public Output {
  public:
    static const std::string keyword;
    
  private:
    string input_filename, output_filename;
    ifstream input_file;
    string output_dir;    
    vector<ofstream*> output_file;      // vector of ofstream is not allowed
    vector<double> output_time;

    double average_time;
    unsigned int output_freq;
    unsigned int time_count;
    unsigned int atom_count;
    unsigned int step_count_per_average_time;
    vector<unsigned int> output_counter;

  public:
    OutputIonTemperatureLS();
    ~OutputIonTemperatureLS() {};
    OutputIonTemperatureLS(const std::string &output_filename, const std::string &input_filename, double average_time, unsigned int output_freq);


  private:
    Output *doMake(const std::vector<Value> &values) const;
    void doInitialize();
    void doRun(int step);
    void doFinalize(int step);
  
  public:
    std::string getIdNoAlias() const {return keyword;}
    void getParameters(std::vector<Parameter> &parameter) const;
  };

}


#endif
