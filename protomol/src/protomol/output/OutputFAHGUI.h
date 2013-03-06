/*  -*- c++ -*-  */
#ifndef PROTOMOL_OUTPUT_FAH_GUI_H
#define PROTOMOL_OUTPUT_FAH_GUI_H

// Updated for standalone GUI
#if defined (HAVE_GUI) || defined (HAVE_LIBFAH)

#include <protomol/output/Output.h>
#include <protomol/base/Timer.h>
#include <string>

#ifdef HAVE_LIBFAH
namespace FAH {
  class GUIServer;
}
#else
#include <protomol/output/GUIServer.h>
#endif

namespace ProtoMol {
  class Configuration;

  class OutputFAHGUI : public Output {
  public:
    static const std::string keyword;

  protected:
    std::string name;
    int port, portRange;
    std::string projName;
    Real timeout;
    bool pause;
    Timer guiTimer;

#ifdef HAVE_LIBFAH
    FAH::GUIServer *server;
#else
    GUIServer *server;
#endif

  public:
    OutputFAHGUI();
    OutputFAHGUI(const std::string &name, int freq, int port,
                 int prange, const std::string &projn, Real timeout,
                 bool pause);

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
    bool adjustWithDefaultParameters(std::vector<Value> &values,
                                     const Configuration *config) const;

  private:
    void setCoords();
    void setBonds();
    void setAtoms();
  };
}

#endif //  HAVE_GUI
#endif //  PROTOMOL_OUTPUT_FAH_GUI_H
