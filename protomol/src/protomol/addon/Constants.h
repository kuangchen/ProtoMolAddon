#ifndef __MY_CONSTANT_H
#define __MY_CONSTANT_H

#include <protomol/base/PMConstants.h>

using namespace ProtoMol::Constant;

namespace ProtoMolAddon {
  namespace Constant {

    namespace ToSI {
      const double velocity = 1e-10 * SI::TIME_FS / TIMEFACTOR;
      const double position = 1.0 / SI::LENGTH_AA;
      const double time = 1.0 / SI::TIME_FS;
      const double force = SI::KCAL * SI::AVOGADRO * 1e-10;
      const double mass = SI::AMU;
      const double charge = SI::ELECTRON_CHARGE / SQRTCOULOMBCONSTANT;
    }

    const double epsilon_0 = 8.85e-12;
    
    const double force_conv = SI::KCAL * SI::AVOGADRO * 1e-10;
    const double position_conv = 1.0 / SI::LENGTH_AA;

    const double FORCE_CONV = SI::KCAL * SI::AVOGADRO * 1e-10;
    const double POSITION_CONV = 1.0 / SI::LENGTH_AA;
    const double TIME_CONV = 1.0 / SI::TIME_FS;
    const double VELOCITY_CONV = 1e-10 * SI::TIME_FS / TIMEFACTOR;
    const double CHARGE_CONV = 1.0 / SQRTCOULOMBCONSTANT;
    const double HBAR = 1.05e-34;
    const double EPSILON_0 = 8.85e-12;
  }
}

#endif
