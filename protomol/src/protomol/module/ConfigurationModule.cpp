#include <protomol/module/ConfigurationModule.h>
#include <protomol/config/Configuration.h>
#include <protomol/io/ConfigurationReader.h>

#include <protomol/config/CommandLine.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/config/InputValue.h>

#include <iostream>
#include <string>

using namespace std;
using namespace ProtoMol;

defineInputValueAndText(InputConfig, "config", "The configuration file name.");

void ConfigurationModule::init(ProtoMolApp *app) {
  config = &app->config;

  CommandLineOption *option;
  CommandLineOption::ActionBase *action;
  CommandLine &cmdLine = app->cmdLine;

  // Keywords
  action = new CommandLineOption::
             Action<ConfigurationModule>(this,
                                         &ConfigurationModule::listKeywords);
  option = cmdLine.add(0, "keywords", action,
    "List all available keywords and exit.");

  // Config
  action = new CommandLineOption::
             Action<ConfigurationModule>(this, &ConfigurationModule::configure);
  option = cmdLine.add(0, "config", action, "Set the configuration file.");
  option->addRequiredArg("filename");

  // Add system keywords
  InputConfig::registerConfiguration(config);
}

int ConfigurationModule::listKeywords() {
  const Configuration *config = ConfigurationModule::config;

  Configuration::const_iterator it;
  for (it = config->begin(); it != config->end(); it++) {
    string key = it->first;
    string type = it->second.getDefinitionTypeString();
    string help = config->getText(key);

    cout << "  --" << key << " " << type << " ";

    fillFormat(cout, help, 6 + key.length() + type.length(), 40);
  }

  return -1;
}

int ConfigurationModule::configure(const vector<string> &args) {
  string configfile = args[1];

  if (config->valid(InputConfig::keyword))
    THROW(string("Configuration file already set to '") +
      config->get(InputConfig::keyword).getString() + "'.");

  config->set(InputConfig::keyword, configfile);

  // Read configuration file
  ConfigurationReader configReader;
  if (!configReader.open(configfile))
    THROW(string("Can't open configuration file '") + configfile + "'.");

  if (!(configReader >> *config))
    THROW(string("Could not read configuration file '") + configfile + "'.");

  return 0;
}
