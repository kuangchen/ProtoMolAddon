#ifndef __CROSS_SECTION_H
#define __CROSS_SECTION_H

#include <memory>
#include <string>
#include <random>
#include <utility>

namespace ProtoMolAddon {
  namespace SympatheticCooling {

    using namespace std;
    
    class DiscreteCrossSection {
    public:
      vector<double> theta;
      vector<double> sigma;

    private:
      random_device rd;
      shared_ptr< piecewise_constant_distribution<double> > theta_dice;
      shared_ptr< uniform_real_distribution<double> > phi_dice;

    public:
      DiscreteCrossSection(const string &fname);
      DiscreteCrossSection();
      pair<double, double> ResampleSolidAngle() const;
    };

  }
}

#endif
