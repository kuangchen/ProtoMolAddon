#ifndef PROTOMOL_CHECKPOINT_MODULE_H
#define PROTOMOL_CHECKPOINT_MODULE_H

#include <protomol/base/Module.h>
#include <protomol/config/InputValue.h>

#include <string>

namespace ProtoMol {
  class CheckpointModule : public Module {
    bool enabled;

  public:
    CheckpointModule() : enabled(false) {}

    const std::string getName() const {return "Checkpoint";}
    int getPriority() const {return -10;}
    const std::string getHelp() const {return "Provides checkpointing support";}

    void init(ProtoMolApp *app);
    void configure(ProtoMolApp *app);
    void read(ProtoMolApp *app);
    void postBuild(ProtoMolApp *app);

    std::string WithoutExt(const std::string &path);
  };
}

#endif // PROTOMOL_CHECKPOINT_MODULE_H

