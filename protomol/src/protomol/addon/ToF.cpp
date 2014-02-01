#include <protomol/addon/ToF.h>
#include <protomol/addon/LuaConfigReader.h>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <boost/range.hpp>
#include <boost/multi_array.hpp>
#include <vector>
#include <iomanip>

using namespace std;
using namespace ProtoMol;
using namespace ProtoMolAddon;
using namespace ProtoMolAddon::Lua;

ToF::ToF() {}

ToF::ToF(const string &def): elct(12) {
  LuaConfigReader reader(def);

  for (int i=0; i<12; i++) {
    ostringstream oss;
    oss << "tof.el" << i;
    elct[i].SetLabel(oss.str());
    
    try {
      string volt_file = reader.GetValue<string>((oss.str() + ".volt").c_str());
      ifstream v(volt_file);
      v >> elct[i].Volt();
    } catch (LuaConfigReaderException &e) {}
    
    try {
      string potl_file = reader.GetValue<string>((oss.str() + ".potl").c_str());
      ifstream p(potl_file);
      p >> elct[i].Potl();
    } catch (LuaConfigReaderException &e) {}
  }

  for (int i=0; i<12; i++) {
    ostringstream oss;
    oss << "tof.el" << i;

    if (!elct[i].Volt().Initialized()) {
      string mirror_el = reader.GetValue<string>((oss.str() + ".volt_copy_from").c_str());
      int k = stoi(mirror_el.substr(2, 3));
      elct[i].Volt(elct[k].Volt());
    }

    if (!elct[i].Potl().Initialized()) {
      string mirror_el2 = reader.GetValue<string>((oss.str() + ".potl_mirror_from").c_str());
      string reflection = reader.GetValue<string>((oss.str() + ".potl_reflection").c_str());
      int k = stoi(mirror_el2.substr(2, 3));
      elct[i].Potl(elct[k].Potl()).SetReflection(reflection);
    }
  } 
}

namespace ProtoMolAddon {
  ostream& operator<<(ostream& os, const ToF& tof) {
    for (auto &el: tof.elct) os << el;
    return os;
  }
}



double ToF::GetTotalRealTimePotential(const Vector3D& pos, double t, const boost::array<int, 3>& offset=boost::array<int, 3>()) {
  vector<double> potl;
  transform(elct.begin(), elct.end(), potl.begin(), 
	    [&pos, t, &offset](const Electrode& r) { return r.GetRealTimePotential(pos, t, offset); } );

  return accumulate(potl.begin(), potl.end(), 0.0);

  // double total_potl = 0;
  // //cout << "pos = " << pos << "\n";
  // for (auto& el: elct) {
  //   double a = el.GetRealTimePotential(pos, t, offset);
  //   //cout << el.GetLabel() << "\t" << a << "\n";
  //   total_potl += a;
  // }

  // return total_potl;
}

double ToF::GetTotalRealTimeInterpolatedPotential(const Vector3D& pos, double t) {
  vector<double> potl;
  transform(elct.begin(), elct.end(), potl.begin(), 
	    [&pos, t](const Electrode& r) { return r.GetRealTimeInterpolatedPotential(pos, t); } );

  return accumulate(potl.begin(), potl.end(), 0.0);

}

void ToF::GetForce(double charge, const Vector3D &pos, double t, Vector3D& force) {
  Vector3D gradient;

  for (int i=0; i<3; i++) {
    boost::array<int, 3> offset_plus = {0,0,0};
    offset_plus[i] = 1; 
    double total_pot_plus = GetTotalRealTimePotential(pos, t, offset_plus);
    double total_pot_minus = GetTotalRealTimePotential(pos, t);
    
    gradient[i] = -(total_pot_plus - total_pot_minus)/(elct[0].Potl().GetDx()[i]);
  }
  force = gradient * charge;
}  

ToF::~ToF() {}

void ToF::Test() {
  ToF tof("test.lua");
  Vector3D pos_raw[] = { Vector3D(0, 0, 0), Vector3D(5e-3, 0, 0), Vector3D(10e-3, 0, 0), 
			 Vector3D(15e-3, 0, 0), Vector3D(20e-3, 0, 0), Vector3D(0, 3e-3, 0), 
			 Vector3D(0, -5e-3, 0), Vector3D(0, 0, 5e-3), Vector3D(0, 0, 15e-3),
			 Vector3D(0, 0, -8e-3), Vector3D(-4e-3, 0, 0), Vector3D(1e-3, 2e-3, 3e-3),
			 Vector3D(3e-3, 3e-3, 3e-3), Vector3D(1, 0, 0) };
  vector<Vector3D> pos(pos_raw, pos_raw+14);

  double real_potl_raw[] = {0.287105, 0.286050, 0.274284, 0.228074, 
			   0.167541, 0.342868, 0.206029, 0.270524, 
			   0.125447, 0.284812, 0.290951, 0.312909, 
			    0.318143, 0};
  vector<double> real_potl(real_potl_raw, real_potl_raw+14);

  vector<double> sim_potl;
  transform(pos.begin(), pos.end(), back_inserter(sim_potl), 
	    [&tof](const Vector3D& p) { return tof.GetTotalRealTimeInterpolatedPotential(p, 0);} );
                                                  
  cout << fixed << setprecision(4);

  for (int i=0; i<14; i++)
    cout << real_potl[i] << "\t" << sim_potl[i] << "\t" << fabs(sim_potl[i]/real_potl[i]-1) << "\n";
  
  Vector3D origin(0, 0, 1e-5);
  Vector3D gradient;

  for (int i=0; i<3; i++) {
    boost::array<int, 3> offset_plus = {0,0,0};
    boost::array<int, 3> offset_minus = {0,0,0};

    offset_plus[i] = 1; 

    double total_pot_plus = 0;
    double total_pot_minus = 0;
    for (auto &el: tof.elct) {
      total_pot_plus += el.Potl().GetPotential(origin, offset_plus);
      total_pot_minus += el.Potl().GetPotential(origin, offset_minus);
    }
    gradient[i] = -(total_pot_plus - total_pot_minus)/(tof.elct[0].Potl().GetDx()[i]);
  }
  cout << "origin = " << origin << "\n";
  cout << "gradient = " << gradient << "\n";
}
