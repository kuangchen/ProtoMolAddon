#include <protomol/addon/sympathetic_cooling/CrossSection.h>
#include <cmath>
#include <fstream>
#include <stdexcept>

using namespace ProtoMolAddon::SympatheticCooling;

random_device CrossSection::rd;

CrossSection::CrossSection(const string &fname) 
{
  ifstream is(fname.c_str());
  
  if (!is)
    throw runtime_error(string("Error opening cross-section file ") + fname);

  double t = 0, s = 0, t_prev = 0;
  vector<double> theta_interp;

  // Add exception for bad data format?
  while (is >> t >> s) {
    if (!theta.empty()) {
      theta.push_back((t+t_prev)/2);
      sigma.push_back(s*(t-t_prev)/2);
    }
    t_prev = t;
  }

  theta_dice.reset(new piecewise_constant_distribution<double>(theta.begin(), 
							       theta.end(),
							       sigma.begin()+1));

  phi_dice.reset(new uniform_real_distribution<double>(0, 2*M_PI));
}
  
pair<double, double> CrossSection::ResampleSolidAngle() const {
  return make_pair<double, double>((*theta_dice)(rd), (*phi_dice)(rd));
}

CrossSection::CrossSection() : 
  theta_dice(),
  phi_dice()  
{

}

//CrossSection::
