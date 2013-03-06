#include <protomol/io/CheckpointConfigReader.h>

#include <protomol/base/Report.h>
#include <protomol/base/StringUtilities.h>

#include <protomol/io/DCDTrajectoryReader.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;


CheckpointConfigReader::CheckpointConfigReader() : Reader() {}


CheckpointConfigReader::CheckpointConfigReader(const string &filename) :
  Reader(filename) {}


bool CheckpointConfigReader::tryFormat() {
  if (!open()) return false;

  string header = getline();
    
  cout << header << endl;

  return header == "!Checkpoint File!";
}


bool CheckpointConfigReader::readBase(Configuration &conf, Random &rand) {
  if (!tryFormat()) {
    cout << "Invalid checkpoint" << endl;
    return false;
  }
    
  cout << "Reading checkpoint" << endl;

  int id = 0, step = 0;
    
  string line;
  while (std::getline(file, line)) {
    if (line.find("#ID") != string::npos) file >> id;
    if (line.find("#Step") != string::npos) file >> step;
    if (line.find("#Random") != string::npos) file >> rand;
  }

  // Update initial checkpoint perameters
  conf["CheckpointStart"] = id + 1;

  // Update position file
  conf["posfile"] = Append(conf["CheckpointPosBase"], id) + ".pos";

  // Update velocities file
  conf["velfile"] = Append(conf["CheckpointVelBase"], id) + ".vel";

  // Update energy file
  if (conf.valid("allEnergiesFile"))
    conf["allEnergiesFile"] = Append(conf["allEnergiesFile"], id);

  // Update firststep
  int firststep = toInt(conf["firststep"]);
  conf["firststep"] = step;

  // Update total steps
  int numsteps = toInt(conf["numsteps"]) - (step - firststep);
  if (numsteps < 0) numsteps = 0;
  conf["numsteps"] = numsteps;

  // DCD file in use?
  if( conf.valid("DCDFile") ){

      //test for file exists and recover 'firstframe'
      const std::string temp = conf["DCDFile"];
      DCDTrajectoryReader dcdTestRead( temp );

      //sucessfully opened?
      //If not then Protomol will just pen a new DCD with the updated 'firststep'
      if( dcdTestRead.tryFormat() ){

          int dcdFirststep = dcdTestRead.readFirstStep();

          report << plain << "'DCD to append to' first step was " << dcdFirststep << "." << endr;

          //valid first step?
          //if not start again!
          if( dcdFirststep >= 0 ){
          
              //get DCD write frequency
              int dcdFrequency = 1; //if not set default is 1

              //set explicitly?
              if( conf.valid("DCDFileOutputFreq")) {
                  dcdFrequency = toInt(conf["DCDFileOutputFreq"]);
              }else{    //else use default output frequency
                  if(conf.valid("outputfreq")){
                    dcdFrequency = toInt(conf["outputfreq"]);
                  }
              }

              //set frame offset value, number of frames stored since simulation began
              //dcdFrequency guarenteed +ve by parser
              //add 1 for first frame output
              conf["DCDFileFrameOffset"] = (step - dcdFirststep) / dcdFrequency + 1;
                                            //(step - firststep) / dcdFrequency + 1;

          }

      }
      
  }

  return !file.fail();
}


bool CheckpointConfigReader::readIntegrator(Integrator *integ) {
  if (!tryFormat()) return false;
    
  string line;
  while (std::getline(file, line))
    if (line.find("#Integrator") != string::npos)
      file >> *integ;

  return !file.fail();
}
