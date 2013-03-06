/*  -*- c++ -*-  */
#ifndef CHECKPOINTCONFIGWRITER_H
#define CHECKPOINTCONFIGWRITER_H

#include <protomol/io/Writer.h>
#include <protomol/base/Random.h>
#include <protomol/integrator/Integrator.h>

namespace ProtoMol {
  //____CheckpointConfigWriter

  class CheckpointConfigWriter : public Writer {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    CheckpointConfigWriter();
    explicit CheckpointConfigWriter(const std::string &filename);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CheckpointConfig
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    bool write(const int& id, const int& steps, const Random& rand,
               const Integrator* integ);
  };
}

#endif /* CheckpointConfigWRITER_H */
