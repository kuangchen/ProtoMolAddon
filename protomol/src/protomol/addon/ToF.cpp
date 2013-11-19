#include <protomol/addon/ToF.h>
#include <protomol/addon/RealElectrode.h>
#include <protomol/addon/MirrorElectrode.h>
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

ToF::ToF(const string &def) {
  electrode.resize(4);
  for (auto& el: electrode) 
    el.resize(3);

  LuaConfigReader reader(def);
  
  for (int i=0; i<4; i++)
    for (int j=0; j<3; j++) {
      
      ostringstream oss;
      oss << "tof.el" << i+1 << j+1;

      try {
	//cout << "i = " << i <<"\tj = " << j << "\n";
	string volt_file = reader.GetValue<string>((oss.str() + ".volt").c_str());
	string pot_file = reader.GetValue<string>((oss.str() + ".pot").c_str());

  	ostringstream label_stream;
  	electrode[i][j] = new RealElectrode(label_stream.str());
            
  	(*static_cast<RealElectrode *>(electrode[i][j])).AssignVoltageFromFile(volt_file).AssignPotentialFromFile(pot_file);

      } catch (LuaConfigReaderException &e) {
      }
    }


  for (int i=0; i<4; i++)
    for (int j=0; j<3; j++) {
      if (electrode[i][j]) continue;
      ostringstream oss;
      oss << "tof.el" << i+1 << j+1;

      string mirror_el = reader.GetValue<string>((oss.str() + ".mirror").c_str());
      string reflection = reader.GetValue<string>((oss.str() + ".reflection").c_str());
      
      int k = std::stoi(mirror_el.substr(2, 1))-1;
      int l = std::stoi(mirror_el.substr(3, 1))-1;

      ostringstream label_stream;
      //      cout << mirror_el << " k = " << k << "\tl = " << l <<"\n";
      assert(electrode[k][l]);
	
      electrode[i][j] = new MirrorElectrode(label_stream.str(), static_cast<RealElectrode*>(electrode[k][l]), reflection);
      
    }

}

void ToF::DumpElectrodeInfo(ostream& os) {
  for (int i=0; i<4; i++)
    for (int j=0; j<3; j++) {
      os << i << "\t" <<j << "\n";
      if (electrode[i][j])
	  electrode[i][j] -> DumpInfo(os);
    
      else 
	os << i << "\t" <<j << "\t empty" << "\n";
    }
}

double ToF::GetTotalPotential(const Vector3D& pos, double t, const array<int, 3>& offset) {
  double total_pot = 0;

  for (auto& rod: electrode)
    for (auto& seg: rod) 
      total_pot += seg->GetNNVoltage(t, 0) * seg->GetNNPotential(pos, offset);

  return total_pot;
}

void ToF::GetForce(double charge, const Vector3D &pos, double t, Vector3D& force) {
  for (int i=0; i<3; i++) {
    std::array<int, 3> offset_plus = {0,0,0}, offset_minus = {0,0,0};
    offset_plus[i] = 1; 
    offset_minus[i] = -1;
    double total_pot_plus = GetTotalPotential(pos, t, offset_plus);
    double total_pot_minus = GetTotalPotential(pos, t, offset_minus);
    
    force[i] = (total_pot_plus - total_pot_minus)/electrode[0][0]->GetDx()[i];
  }

  force *= charge;
}

ToF::~ToF() {
  for (auto& rod: electrode)
    for (auto& seg: rod) {
      if (seg) delete seg;
    }
}

void ToF::Test() {
  ToF tof("test.lua");
  vector<Vector3D> pos;

  pos.push_back(Vector3D(0, 0, 0));
  pos.push_back(Vector3D(5e-3, 0, 0));
  pos.push_back(Vector3D(10e-3, 0, 0));
  pos.push_back(Vector3D(15e-3, 0, 0));
  pos.push_back(Vector3D(20e-3, 0, 0));
  pos.push_back(Vector3D(0, 3e-3, 0));
  pos.push_back(Vector3D(0, -5e-3, 0));
  pos.push_back(Vector3D(0, 0, 5e-3));
  pos.push_back(Vector3D(0, 0, 15e-3));
  pos.push_back(Vector3D(0, 0, -8e-3));
  pos.push_back(Vector3D(-4e-3, 0, 0));
  pos.push_back(Vector3D(1e-3, 2e-3, 3e-3));
  pos.push_back(Vector3D(3e-3, 3e-3, 3e-3));

  double voltage[4][3] = { {0.667, 0.333, 0.44}, 
			   {-0.5, 0.25, 0.125},
			   {0.2, 0.1, 0.55},
			   {-0.2, 0.5, -0.1} };
  
  vector<double> real_pot;
  real_pot.push_back(0.287105);
  real_pot.push_back(0.286050);
  real_pot.push_back(0.274284);
  real_pot.push_back(0.228074);
  real_pot.push_back(0.167541);
  real_pot.push_back(0.342868);
  real_pot.push_back(0.206029);
  real_pot.push_back(0.270524);
  real_pot.push_back(0.125447);
  real_pot.push_back(0.284812);
  real_pot.push_back(0.290951);
  real_pot.push_back(0.312909);
  real_pot.push_back(0.318143);

  vector<double> sim_pot;
  vector<double> variation;
  for (auto& p: pos) {
    double total_pot = 0;
    double var;
    for (int i=0; i<4; i++)
      for (int j=0; j<3; j++) {
	array<double, 3> f;
	tof.electrode[i][j]->GetFraction(p, f);
	double ave_pot = tof.electrode[i][j]->GetNNPotential(p, {0, 0, 0}) * (1-f[0]) * (1-f[1]) * (1-f[2]) +
	  tof.electrode[i][j]->GetNNPotential(p, {0, 0, 1}) * (1-f[0]) * (1-f[1]) * f[2] +
	  tof.electrode[i][j]->GetNNPotential(p, {0, 1, 0}) * (1-f[0]) * f[1] * (1-f[2]) +
	  tof.electrode[i][j]->GetNNPotential(p, {0, 1, 1}) * (1-f[0]) * f[1] * f[2] +
	  tof.electrode[i][j]->GetNNPotential(p, {1, 0, 0}) * f[0] * (1-f[1]) * (1-f[2]) +
	  tof.electrode[i][j]->GetNNPotential(p, {1, 0, 1}) * f[0] * (1-f[1]) * f[2] +
	  tof.electrode[i][j]->GetNNPotential(p, {1, 1, 0}) * f[0] * f[1] * (1-f[2]) +
	  tof.electrode[i][j]->GetNNPotential(p, {1, 1, 1}) * f[0] * f[1] * f[2];
	
	vector<double> nn;
	nn.push_back(tof.electrode[i][j]->GetNNPotential(p, {0, 0, 0}));
	nn.push_back(tof.electrode[i][j]->GetNNPotential(p, {0, 0, 1}));
	nn.push_back(tof.electrode[i][j]->GetNNPotential(p, {0, 1, 0}));
	nn.push_back(tof.electrode[i][j]->GetNNPotential(p, {0, 1, 1}));
	nn.push_back(tof.electrode[i][j]->GetNNPotential(p, {1, 0, 0}));
	nn.push_back(tof.electrode[i][j]->GetNNPotential(p, {1, 0, 1}));
	nn.push_back(tof.electrode[i][j]->GetNNPotential(p, {1, 1, 0}));
	nn.push_back(tof.electrode[i][j]->GetNNPotential(p, {1, 1, 1}));

	double min_val = *min_element(nn.begin(), nn.end());
	double max_val = *max_element(nn.begin(), nn.end());
	var = max_val / min_val - 1;
	total_pot += ave_pot * voltage[i][j];	
      }
    variation.push_back(var);
    sim_pot.push_back(total_pot);
  }

  cout << fixed << setprecision(4);
  for (int i=0; i<13; i++)
    cout << real_pot[i]  << "\t" << sim_pot[i] << "\t" << fabs(sim_pot[i]/real_pot[i]-1) << "\t" << variation[i] << "\n";
  
  
}
