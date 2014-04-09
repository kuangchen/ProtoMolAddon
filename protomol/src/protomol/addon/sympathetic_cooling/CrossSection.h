#ifndef __CROSS_SECTION_H
#define __CROSS_SECTION_H

#include <memory>
#include <string>
#include <random>
#include <utility>

namespace ProtoMolAddon {
  namespace SympatheticCooling {

    using namespace std;
    
    class CrossSection {
    public:
      vector<double> theta;
      vector<double> sigma;

    private:
      static random_device rd;
      shared_ptr< piecewise_constant_distribution<double> > theta_dice;
      shared_ptr< uniform_real_distribution<double> > phi_dice;

    public:
      CrossSection(const string &fname);
      CrossSection();
      pair<double, double> ResampleSolidAngle() const;
    };

  }
}

#endif
