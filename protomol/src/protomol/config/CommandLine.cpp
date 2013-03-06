#include <protomol/config/CommandLine.h>

#include <protomol/base/Exception.h>
#include <protomol/type/String.h>
#include <protomol/config/Configuration.h>
#include <protomol/module/ConfigurationModule.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/base/SystemUtilities.h>
#include <protomol/ProtoMolApp.h>

#include <stdlib.h>

using namespace std;
using namespace ProtoMol;

CommandLine::CommandLine(Configuration *config) :
  config(config), name("ProtoMol") {}

void CommandLine::add(CommandLineOption *option) {
  if (!option) THROW("option is NULL");

  CommandLineOption::names_t::iterator it;
  for (it = option->names.begin(); it != option->names.end(); it++) {
    if (optionMap.find(*it) != optionMap.end())
      THROW("option '" + *it + "' already exists");

    optionMap[*it] = option;
  }

  options.push_back(option);
}

CommandLineOption *CommandLine::add(const char shortName, const string longName,
                                    CommandLineOption::ActionBase *action,
                                    const string help) {
  CommandLineOption *option =
    new CommandLineOption(shortName, longName, action, help);
  add(option);

  return option;
}

int CommandLine::parse(int argc, char *argv[]) {
  return parse(vector<string>(argv, argv + argc));
}

int CommandLine::parse(const vector<string> &args) {
  name = args[0];
  bool configSet = false;

  for (unsigned int i = 1; i < args.size();) {
    string optionStr = args[i];
    CommandLineOption *option = optionMap[optionStr];

    if (option) {
      vector<string> cmdArgs;

      // Self arg
      cmdArgs.push_back(args[i++]);

      // Required args
      for (unsigned int j = 0; j < option->requiredArgs.size(); j++) {
        if (args.size() <= i)
          THROW(string("Missing required argument for option ") + optionStr);

        cmdArgs.push_back(args[i++]);
      }

      // Optional args
      for (unsigned int j = 0; j < option->optionalArgs.size(); j++) {
        if (args.size() <= i || args[i][0] == '-') break;

        cmdArgs.push_back(args[i++]);
      }

      // Call action
      if (option->action) {
        int ret = (*option->action)(cmdArgs);
        if (ret == -1) return -1;
      }

    } else if (args[i][0] == '-' && args[i][1] == '-') {
      string key = &args[i++][2];

      const Configuration *const_config = config;
      Configuration::const_iterator it = const_config->find(key);
      if (it == const_config->end())
        THROW(string("Invalid keyword '") + key + "'");

      if (it->second.getType() == ValueType::VECTOR ||
          it->second.getType() == ValueType::VECTOR3D)
        THROW(string("Keyword '") + key +
          "' with Vector type cannot be set from the command line.");

      if (i == args.size())
        THROW(string("Missing argument for keyword '") + key + "'.");
      string val = args[i++];

      if (!config->set(key, val))
        THROW(string("Invalid value '") + val + "' for keyword '" + key + "'.");

    } else if (!configSet) {
      vector<string> cmdArgs;

      cmdArgs.push_back("--config");
      cmdArgs.push_back(args[i++]);

      (*optionMap["--config"]->action)(cmdArgs);
      configSet = true;
    } else THROW(string("Invalid argument '") + args[i] + "'");
  }

  return 0;
}

int CommandLine::usageAction(const vector<string> &args) {
  usage(cerr, name);
  return -1;
}

void CommandLine::usage(ostream &stream, const string &name) {
  stream
    << "Usage: " << name << " [--config] <filename> [--option args]..." << endl
    << "Options:" << endl;

  for (unsigned int i = 0; i < options.size(); i++) {
    int count = 2;
    stream << "  ";

    // Short option names
    bool first = true;
    for (CommandLineOption::names_t::iterator it = options[i]->names.begin();
         it != options[i]->names.end(); it++) {
      if ((*it)[0] == '-' && (*it)[1] == '-') continue;

      if (first) first = false;
      else {
        stream << "|";
        count++;
      }

      stream << *it;
      count += (*it).length();
    }

    // Long option names
    for (CommandLineOption::names_t::iterator it = options[i]->names.begin();
         it != options[i]->names.end(); it++) {
      if ((*it)[0] == '-' && (*it)[1] != '-') continue;

      if (first) first = false;
      else {
        stream << "|";
        count++;
      }

      stream << *it;
      count += (*it).length();
    }

    // Args
    for (unsigned int j = 0; j < options[i]->requiredArgs.size(); j++) {
      stream << " <" << options[i]->requiredArgs[j] << ">";
      count += 3 + options[i]->requiredArgs[j].length();
    }

    for (unsigned int j = 0; j < options[i]->optionalArgs.size(); j++) {
      stream << " [" << options[i]->optionalArgs[j] << "]";
      count += 3 + options[i]->optionalArgs[j].length();
    }

    // Help
    if (count >= 40) {
      count = 0;
      stream << endl;
    }

    fillFormat(stream, options[i]->help, count, 40);
  }
}

int CommandLine::splashAction(const std::vector<std::string> &args) {
  ProtoMolApp::splash(cout);
  return -1;
}

#ifdef DEBUG
int CommandLine::enableStackTraceAction(const vector<string> &args) {
#ifdef HAVE_STACK_TRACE
  Debugger::initStackTrace(getCanonicalPath(name));
#endif

  return 0;
}

#endif
