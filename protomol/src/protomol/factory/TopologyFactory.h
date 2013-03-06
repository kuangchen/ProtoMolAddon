/* -*- c++ -*- */
#ifndef TOPOLOGYFACTORYDETAILS_H
#define TOPOLOGYFACTORYDETAILS_H

#include <protomol/base/Factory.h>
#include <protomol/topology/GenericTopology.h>

namespace ProtoMol {
  class Configuration;

  //________________________________________ TopologyFactory
  class TopologyFactory : public Factory<GenericTopology> {
  public:
    TopologyFactory() {}
    virtual ~TopologyFactory() {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From Factory<GenericTopology>
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    virtual void registerHelpText() const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class TopologyFactory
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void registerAllExemplarsConfiguration(Configuration *config) const;
    GenericTopology *make(const Configuration *config) const;
    GenericTopology *make(const std::string &id,
                          const std::vector<Value> &values) const;
  };
}
#endif /* TOPOLOGYFACTORYDETAILS_H */
