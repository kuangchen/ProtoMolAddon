#ifndef __CROSS_SECTION_H_
#define __CROSS_SECTION_H_

#include <iostream>
#include <memory>
#include <string>
#include <random>
#include <utility>
#include <vector>

namespace ProtoMolAddon {
  namespace SympatheticCooling {

    using namespace std;
    
    class DiscreteCrossSection {
    private:
      struct CrossSectionSpec {

	struct ThetaSigmaPair {
	  double theta;
	  double sigma;
	  
	  friend istream& operator>> (istream &is, ThetaSigmaPair &ts) { is >> ts.theta >> ts.sigma; return is; }
	  friend bool operator< (const ThetaSigmaPair &p1, const ThetaSigmaPair &p2) { return p1.theta < p2.theta; }
	};
	
	vector<ThetaSigmaPair> theta_sigma_array;
	
	CrossSectionSpec();
	CrossSectionSpec(const std::string &fname);
      };

    public:
      CrossSectionSpec spec;
      random_device rd;
      shared_ptr< piecewise_constant_distribution<double> > theta_dice;
      shared_ptr< uniform_real_distribution<double> > phi_dice;

      DiscreteCrossSection();
      
      DiscreteCrossSection(const string &fname);
      pair<double, double> ResampleSolidAngle() const;
      
      friend bool operator< (const DiscreteCrossSection &cs1, const DiscreteCrossSection &cs2);
    };

  }
}

#endif
