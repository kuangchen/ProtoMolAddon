#include <protomol/addon/sympathetic_cooling/DiscreteCrossSection.h>
#include <boost/program_options.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iterator>
#include <stdexcept>

namespace po = boost::program_options;
using namespace ProtoMolAddon::SympatheticCooling;

DiscreteCrossSection::CrossSectionSpec::CrossSectionSpec() : theta_sigma_array(0) {
}

DiscreteCrossSection::CrossSectionSpec::CrossSectionSpec(const std::string &fname) 
  : theta_sigma_array(0)
{
  typedef DiscreteCrossSection::CrossSectionSpec::ThetaSigmaPair ts_pair;
  
  ifstream is(fname.c_str());
  
  if (!is)
    throw runtime_error(string("Error opening cross-section file ") + fname);

  po::variables_map vm; 
  po::options_description desc("CrossSection Spec"); 
  desc.add_options()
    ("CrossSection.Value",  po::value< std::vector<ts_pair> >(&theta_sigma_array)->required(), "Cross-section Value");

  po::store(po::parse_config_file(is, desc, true), vm); 
  po::notify(vm);
  
  sort(theta_sigma_array.begin(), theta_sigma_array.end());
  
  if (theta_sigma_array.front().theta < 0 || theta_sigma_array.back().theta > M_PI) 
    throw runtime_error("Invalid input theta range");
}  
  
DiscreteCrossSection::DiscreteCrossSection() {}

DiscreteCrossSection::DiscreteCrossSection(const std::string &fname) 
: spec(fname), 
  theta_dice(NULL), 
  phi_dice(new uniform_real_distribution<double>(0, 2*M_PI)) 
{
  
  std::vector<double> theta;
  std::vector<double> weight;

  double prev = 0;  
  double curr;

  for (unsigned int i=0; i<spec.theta_sigma_array.size(); i++) {
    curr = spec.theta_sigma_array[i].theta;

    double middle = (curr + prev) / 2.0;    
    theta.push_back(middle);

    weight.push_back(spec.theta_sigma_array[i].sigma * sin(middle) * (curr-prev));
    prev = curr;
  }

  theta.push_back((prev + M_PI)/2.0);
 
  theta_dice.reset(new piecewise_constant_distribution<double>(theta.begin(),
							       theta.end(),
							       weight.begin()));
}

pair<double, double> DiscreteCrossSection::ResampleSolidAngle() const {
  return make_pair<double, double>((*theta_dice)(rd), (*phi_dice)(rd));
}

