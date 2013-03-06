/*  -*- c++ -*-  */
#ifndef PROTOMOL_OUTPUT_H
#define PROTOMOL_OUTPUT_H

#include <protomol/base/Makeable.h>

namespace ProtoMol {
  class ProtoMolApp;

  /***
     Base class of all Output classes to dump data at a given
     frequency.  The actual output frequency is defined by global
     output frequency times the the frequency of the concret
     class. The global output frequency is also  used to define th
     number of steps to run the integrator.  It keeps pointers to th
     topology, positions, velocities, and integrator.  Output object
     are aggregated in OutputCollection, which  invokes them one by
     one at application level.  If you only need to print some output
     to a file you should rather inherit from OutputFile, then inherit
     directly from Output.
   */
  class Output : public Makeable<Output> {
  public:
    static const std::string scope;
    const ProtoMolApp *app;

  protected:
    int firstStep;
    int nextStep;
    int lastStep;
    int outputFreq; // /< Output freqeuncy

  public:
    Output(int freq = 0);

    // / To initialize the object, before the simulation starts.
    virtual void initialize(const ProtoMolApp *app);

    // / Called at each step (e.g., printing total energy on the screen),
    // / takes care of the output frequency.  Returns true if it ran.
    virtual bool run(int step);

    // / At the end of the simulation (e.g., writing final positions), and
    // / calls first run() to ensure that run is called for the last
    // / step, if needed.
    virtual void finalize(int step);

    // / Factory method to create a complete output object from its prototy
    virtual Output *make(const std::vector<Value> &values) const;

    // / Should return true if the concrete object is defined/specified in
    // / Configuration by the user. Normally if gedId() has a valid valu
    // / in Configuration.
    virtual bool isIdDefined(const Configuration *config) const;

    // / Defines if the output object supports do<getId()> to enable or disabl
    // / the output.
    virtual bool addDoKeyword() const {return true;}

    int getFirstStep() const {return firstStep;}
    int getLastStep() const {return lastStep;}
    int getOutputFreq() const {return outputFreq;}
    int getNext() const {return nextStep;}

  private:
    // / Hook method of initialize, implemented in the concrete cl
    virtual void doInitialize() = 0;

    // / Hook method of run, implemented in the concrete cl
    virtual void doRun(int step) = 0;

    virtual void doFinalize(int step) {};

    //  From class Makabl
  public:
    virtual std::string getScope() const {return scope;}
    void getParameters(std::vector<Parameter> &parameter) const;
    bool adjustWithDefaultParameters(std::vector<Value> &values,
                                     const Configuration *config) const;
  };
}
#endif //  PROTOMOL_OUTPUT_H
