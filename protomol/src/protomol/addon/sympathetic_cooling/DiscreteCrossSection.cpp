#include <protomol/addon/sympathetic_cooling/DiscreteCrossSection.h>
#include <boost/program_options.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iterator>
#include <stdexcept>

namespace po = boost::program_options;

using namespace std;
using namespace ProtoMolAddon::SympatheticCooling;

DiscreteCrossSection::DiscreteCrossSection() :
  rd(), theta_dice(NULL), phi_dice(NULL)
{}

DiscreteCrossSection::DiscreteCrossSection(const string &fname) :
  energy(0), rd(), theta_dice(NULL), phi_dice(NULL)
{
  ifstream is(fname.c_str());
  
  if (!is)
    throw runtime_error(string("Error opening cross-section file ") + fname);

  po::variables_map vm; 
  po::options_description desc("CrossSection Spec"); 
  desc.add_options()
    ("CrossSection.Energy", po::value<double>(&energy)->required(), "Cross-section Energy")
    ("CrossSection.Value",  po::value<vector<theta_sigma> >(&theta_sigma_array)->required(), "Cross-section Value");

  po::store(po::parse_config_file(is, desc, true), vm); 
  po::notify(vm);
  
  sort(theta_sigma_array.begin(), theta_sigma_array.end());
  
  if (theta_sigma_array.front().theta < 0 || theta_sigma_array.back().theta > M_PI) 
    throw runtime_error("Invalid input theta range");
  
  vector<double> theta;
  vector<double> sigma;

  double theta_prev = 0;  
  double theta_current;
 
  for (unsigned int i=0; i<theta_sigma_array.size(); i++) {
    theta_current = theta_sigma_array[i].theta;
    
    theta.push_back((theta_current + theta_prev)/2);
    theta_prev = theta_current;

    sigma.push_back(theta_sigma_array[i].sigma);
  }

  theta.push_back((theta_prev + M_PI)/2);
  
  theta_dice.reset(new piecewise_constant_distribution<double>(theta.begin(), 
							       theta.end(),
							       sigma.begin()));

  phi_dice.reset(new uniform_real_distribution<double>(0, 2*M_PI));
}


  
pair<double, double> DiscreteCrossSection::ResampleSolidAngle() const {
  return make_pair<double, double>((*theta_dice)(rd), (*phi_dice)(rd));
}

namespace ProtoMolAddon {
  namespace SympatheticCooling {

    bool operator< (const DiscreteCrossSection &cs1, const DiscreteCrossSection &cs2) {
      return cs1.energy < cs2.energy;
    }
  }
}
