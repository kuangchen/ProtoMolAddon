#ifndef MODIFIER_MODULE_H
#define MODIFIER_MODULE_H

#include <protomol/base/Module.h>
#include <protomol/config/InputValue.h>

#include <string>

namespace ProtoMol {
  class ProtoMolApp;

  declareInputValue(InputRemoveAngularMomentum, INT, NOCONSTRAINTS);
  declareInputValue(InputRemoveLinearMomentum, INT, NOCONSTRAINTS);
  declareInputValue(InputRattle, BOOL, NOCONSTRAINTS);
  declareInputValue(InputRattleEpsilon, REAL, NOTNEGATIVE);
  declareInputValue(InputRattleMaxIter, INT, NOTNEGATIVE);
  declareInputValue(InputRattleAll, BOOL, NOCONSTRAINTS);

  declareInputValue(InputShake, BOOL, NOCONSTRAINTS);
  declareInputValue(InputShakeEpsilon, REAL, NOTNEGATIVE);
  declareInputValue(InputShakeMaxIter, INT, NOTNEGATIVE);
  declareInputValue(InputShakeAll, BOOL, NOCONSTRAINTS);

  declareInputValue(InputShadow, BOOL, NOCONSTRAINTS);
  declareInputValue(InputShadowOrder, INT, NOTNEGATIVE);
  declareInputValue(InputShadowFreq, INT, NOTNEGATIVE);

  class ModifierModule : public Module {
  public:
    const std::string getName() const {return "Modifier";}
    int getPriority() const {return 0;}
    const std::string getHelp() const {return "";}

    void init(ProtoMolApp *app);
    void postBuild(ProtoMolApp *app);
    void addModifiers(ProtoMolApp *app);
  };
};

#endif // MODIFIER_MODULE_H
