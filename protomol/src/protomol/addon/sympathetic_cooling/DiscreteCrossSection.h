#ifndef __CROSS_SECTION_H
#define __CROSS_SECTION_H

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
    public:
      
      struct theta_sigma {
	double theta;
	double sigma;
	theta_sigma(double t=0, double s=0): theta(t), sigma(s) {}
	
	friend istream& operator>> (istream &is, theta_sigma &ts) { 
	  is >> ts.theta >> ts.sigma; return is; 
	}
	
	friend ostream& operator<< (ostream &os, const theta_sigma &ts) { 
	  os << ts.theta << ts.sigma; return os; 
	}
	
	friend bool operator< (const theta_sigma &ts1, const theta_sigma &ts2) {
	  return ts1.theta < ts2.theta;
	}
      };

    private:
      double energy;
      vector<theta_sigma> theta_sigma_array;
      
    private:
      random_device rd;
      shared_ptr< piecewise_constant_distribution<double> > theta_dice;
      shared_ptr< uniform_real_distribution<double> > phi_dice;

    public:
      DiscreteCrossSection();
      
      DiscreteCrossSection(const string &fname);
      pair<double, double> ResampleSolidAngle() const;
//      vector<double> ReSampleSolidAngleTest() const;
      
      friend bool operator< (const DiscreteCrossSection &cs1, const DiscreteCrossSection &cs2);
    };

  }
}

#endif
