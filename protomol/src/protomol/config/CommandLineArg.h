#ifndef COMMAND_LINE_ARG_H
#define COMMAND_LINE_ARG_H

#include <string>

namespace ProtoMol {
  class CommandLineArg {
  public:
    std::string name;
    bool optional;

    CommandLineArg(const std::string &name, bool optional = false) :
      name(name), optional(optional) {}
  }
};

#endif // COMMAND_LINE_ARG_H

