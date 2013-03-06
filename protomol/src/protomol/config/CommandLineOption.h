#ifndef COMMAND_LINE_OPTION_H
#define COMMAND_LINE_OPTION_H

#include <string>
#include <set>
#include <vector>
#include <ostream>

namespace ProtoMol {
  class CommandLineOption {
  public:

    // Command line action functor declarations
    class ActionBase {
    public:
      virtual ~ActionBase() {}

      virtual int operator()(const std::vector<std::string> &args) = 0;
    };

    template<class T>
    class Action : public ActionBase {
      T *obj;
      typedef int (T::*fpt_t)(const std::vector<std::string> &args);
      fpt_t fpt;

      typedef int (T::*fpt_noargs_t)();
      fpt_noargs_t fpt_noargs;

    public:
      Action(T *obj, fpt_t fpt) : obj(obj), fpt(fpt), fpt_noargs(0) {}
      Action(T *obj, fpt_noargs_t fpt_noargs) :
        obj(obj), fpt(0), fpt_noargs(fpt_noargs) {}
      virtual ~Action() {}

      virtual int operator()(const std::vector<std::string> &args) {
        if (fpt) return (*obj.*fpt)(args);
        else return (*obj.*fpt_noargs)();
      }
    };

    // Member data
    typedef std::set<std::string> names_t;
    names_t names;

    ActionBase *action;

    std::string help;

    typedef std::vector<std::string> args_t;
    args_t requiredArgs;
    args_t optionalArgs;

    // Constructor
    CommandLineOption(const char shortName = 0,
                      const std::string longName = "",
                      ActionBase *action = 0,
                      const std::string help = "");
    
    ~CommandLineOption();

    // Member functions
    void addLongAlias(const std::string &alias) {names.insert("--" + alias);}
    void addShortAlias(const char alias) {
      names.insert(std::string("-") + alias);
    }
    void addRequiredArg(const std::string &arg) {requiredArgs.push_back(arg);}
    void addOptionalArg(const std::string &arg) {optionalArgs.push_back(arg);}
    bool match(const std::string &arg) const;
  };
};


#endif // COMMAND_LINE_OPTION_H
