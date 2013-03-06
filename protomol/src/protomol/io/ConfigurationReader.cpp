#include <protomol/io/ConfigurationReader.h>

#include <protomol/base/Report.h>
#include <protomol/base/StringUtilities.h>

using namespace std; 
using namespace ProtoMol::Report;
using namespace ProtoMol;
//____ConfigurationReader

ConfigurationReader::ConfigurationReader() : Reader() {}

ConfigurationReader::ConfigurationReader(const string &filename) :
  Reader(filename) {}

bool ConfigurationReader::tryFormat() {
  open();
  return !file.fail();
}

bool ConfigurationReader::read(Configuration &config) {
  if (!tryFormat()) return false;

  // Remove comments and reformat
  stringstream all;
  while (!file.eof() && !file.fail()) {
    string line(getline());
    stringstream ss(string(line.begin(), find(line.begin(), line.end(), '#')));
    string str;

    while (ss >> str)
      all << (all.str().empty() ? "" : " ") << str;
  }

  close();
  if (file.fail()) return false;

  // Nothing to do ...
  if (all.str().empty()) {
    report << warning << "Empty configuration file" << endr;
    return true;
  }

  // First get the keyword and then let Value read from istream ...
  string str;
  string bad;
  bool res = true;
  while (all >> str)
    if (!config.empty(str)) {
      if (!bad.empty()) {
        report << recoverable << "Ignoring:" << bad << endr;
        bad = "";
      }
      ios::pos_type start = all.tellg();
      all >> config[str];
      if (!config[str].valid()) {
        ios::pos_type end = all.tellg();
        if (end > start) {
          streamsize len = static_cast<streamsize>(end - start);
          string tmp(len, ' ');
          all.seekg(start);
          all.read(&(tmp[0]), len);
          report << recoverable << "Could not parse \'" <<
          removeBeginEndBlanks(tmp) << "\' for keyword \'" << str
                 << "\', expecting type " <<
          config[str].getDefinitionTypeString() << "." << endr;
        }
      }
    } else {
      res = false;
      bad += " " + str;
    }

  if (!bad.empty())
    report << recoverable << "Ignoring:" << bad << endr;

  return res;
}

namespace ProtoMol {

  ConfigurationReader &operator>>(ConfigurationReader &configReader,
                                  Configuration &config) {
    configReader.read(config);
    return configReader;
  }
}
