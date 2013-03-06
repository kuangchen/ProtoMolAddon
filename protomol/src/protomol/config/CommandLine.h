#ifndef COMMAND_LINE_H
#define COMMAND_LINE_H

#include <protomol/config/CommandLineOption.h>
#include <protomol/base/SmartPointer.h>

#include <string>
#include <map>
#include <vector>
#include <ostream>

namespace ProtoMol {
  class Configuration;

  class CommandLine {
    Configuration *config;

    typedef std::map<const std::string, CommandLineOption *> optionMap_t;
    optionMap_t optionMap;

    typedef std::vector<SmartPointer<CommandLineOption> > options_t;
    options_t options;

    std::string name;

  public:
    CommandLine() {}
    CommandLine(Configuration *config);

    void add(CommandLineOption *option);
    CommandLineOption *add(const char shortName = 0,
                           const std::string longName = "",
                           CommandLineOption::ActionBase *action = 0,
                           const std::string help = "");

    int parse(int argc, char *argv[]);
    int parse(const std::vector<std::string> &args);

    int usageAction(const std::vector<std::string> &args);
    void usage(std::ostream &stream, const std::string &name);
    int splashAction(const std::vector<std::string> &args);

#ifdef DEBUG
    int enableStackTraceAction(const std::vector<std::string> &args);
#endif
  };
};

#endif // COMMAND_LINE_H
