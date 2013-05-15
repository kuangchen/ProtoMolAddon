#include <protomol/addon/NeutralAtom.h>
#include <protomol/addon/NeutralAtomLookupTable.h>
#include <protomol/base/PMConstants.h>
#include <protomol/addon/Constants.h>
#include <protomol/addon/Util.h>
#include <protomol/addon/ProtoMolIon.h>
#include <cmath>

#include <fstream>

using namespace ProtoMolAddon;
using namespace ProtoMolAddon::Util;
using namespace ProtoMol::Constant;
using namespace ProtoMolAddon::Constant;

NeutralAtom::NeutralAtom(double m, double pl, double t, double d, const Vector3D &pos, const Vector3D &vel) :
  mass(m), polarizability(pl), temperature(t), density(d), position(pos), velocity(vel),
  rd(),
  collision_dice(0, 1),
  phi_dice(0, 2*M_PI) {

  sigma = sqrt(SI::BOLTZMANN*temperature/mass);
 
  LoadThetaCache(10000);
  LoadVnCache(10000);
}
  
NeutralAtom::NeutralAtom(LuaConfigReader &reader):
  mass(reader.GetValue<double>("neutral.mass") * SI::AMU),
  polarizability(reader.GetValue<double>("neutral.polarizability") * 1e-30),
  temperature(reader.GetValue<double>("neutral.temperature")),
  density(reader.GetValue<double>("neutral.density")),
  position(Vector3D()), 
  velocity(Vector3D()),
  rd(),
  collision_dice(0, 1),
  phi_dice(0, 2*M_PI) {

  sigma = sqrt(SI::BOLTZMANN*temperature/mass);
  vn_dice = std::normal_distribution<double>(0, sigma);
 
  LoadFLktb();
  LoadThetaCache();
  LoadVnCache();
  LoadCrossSection(reader);
}

bool NeutralAtom::ShouldCollide(const ProtoMolIon &ion, double dt) {
  double rate = GetThermalAverageCollisionRate(ion);
  double dice = collision_dice(rd);

  return (1-exp(-rate*dt))>dice;
}

bool NeutralAtom::ShouldCollide2(const ProtoMolIon &ion, double dt) {
  double rate = GetCollisionRate(ion);
  double dice = collision_dice(rd);

  return (1-exp(-rate*dt))>dice;
}

  
void NeutralAtom::CollideAll(ProtoMolApp *app, double dt) {
  for (unsigned int i=0; i<app->velocities.size(); i++) {
    ProtoMolIon ion(app, i);
    if (ion.charge > 0 && ShouldCollide(ion, dt)) {
      SampleVelocity(ion);
      Collide(ion);
      ion.UpdateProtoMol(app);
    }
  }  
}

void NeutralAtom::CollideAll2(ProtoMolApp *app, double dt) {
  for (unsigned int i=0; i<app->topology->atoms.size(); i++) {
    ProtoMolIon ion(app, i);
    SampleVelocity2(ion);
    if (ion.charge > 0 && ShouldCollide2(ion, dt)) {
      Collide(ion);
      ion.UpdateProtoMol(app);
    }
  }  
}

double NeutralAtom::GetCollisionRate(const ProtoMolIon &ion) {
  double C4 = polarizability * ion.charge * ion.charge / (4*M_PI*EPSILON_0);
  double mu = ion.mass * mass / (ion.mass+mass);
  Vector3D v_rel = velocity-ion.velocity;
  double E = 0.5 * mu * v_rel.normSquared();
  
  return M_PI * pow(C4*C4*mu/HBAR/HBAR/E,1.0/3)*(1+M_PI*M_PI/16)* v_rel.norm() * density;
}

double NeutralAtom::GetThermalAverageCollisionRate(const ProtoMolIon &ion) {
  double C4 = polarizability * ion.charge * ion.charge / (4*M_PI*EPSILON_0);
  double p = M_PI * pow(C4*C4*sigma*2/HBAR/HBAR,1.0/3)*(1+M_PI*M_PI/16)*density;
  double v_norm_reduced = ion.velocity.norm() / sigma;

  auto it = f_lktb.upper_bound(v_norm_reduced);
  if (it==f_lktb.end())
    return p * pow(v_norm_reduced, 1.0/3);

  else {
    auto it2 = it--;
    double f = (v_norm_reduced - it->first) / (it2->first - it->first);
    return (f * it2->second + (1-f) * it->second) * p;
  }
}

void NeutralAtom::SampleVelocity(const ProtoMolIon &ion) {
  Vector3D vi = ion.velocity;
  double vi_norm = vi.norm();
  
  // Populate weight array with analytical expression for vn
  vn_weight.resize(0);
  for (double v: vn_center) 
    vn_weight.push_back(v*exp(-v*v/2/sigma/sigma)*(pow((v+vi_norm),7.0/3)-pow(fabs(v-vi_norm),7.0/3)));
  
  // Sample vn
  pc_dist vn_dice(begin(vn_edge), end(vn_edge), begin(vn_weight)); 
  double vn = vn_dice(rd);

  // Populate weight array with analytical expression for theta
  theta_weight.resize(0);
  for (double t: theta_center) 
    theta_weight.push_back(pow((vn*vn-2*vn*vi_norm*cos(t)+vi_norm*vi_norm), 1.0/6)*sin(t));
  
  // Sample theta
  pc_dist theta_dice(begin(theta_edge), end(theta_edge), begin(theta_weight));
  double theta = theta_dice(rd);

  // Sample phi
  double phi = phi_dice(rd);
  
  // Build and rotate velocity
  Vector3D velocity1 = BuildVector(vn, theta, phi);
  velocity = Rotate(vi, velocity1);
}

void NeutralAtom::SampleVelocity2(const ProtoMolIon &ion) {
  velocity = Vector3D(vn_dice(rd), vn_dice(rd), vn_dice(rd));
}

void NeutralAtom::Collide(ProtoMolIon &ion) {
  double mi = ion.mass;
  Vector3D vi = ion.velocity;

  double mu = mi * mass / (mi+mass);
  double b1 = mi / (mass + mi);
  double b2 = 1-b1;

  Vector3D v_com, v_rel, v_rel_after1, v_rel_after;
  v_com = vi * b1 + velocity * b2;
  v_rel = vi - velocity;
  double energy_in_kelvin = 0.5 * mu * v_rel.normSquared() / SI::BOLTZMANN;

  double v_rel_norm = v_rel.norm();
  pc_dist &dist_use = SelectThetaDist(energy_in_kelvin);

  double theta = dist_use(rd);
  double phi = phi_dice(rd);
  
  v_rel_after1 = BuildVector(v_rel_norm, theta, phi);
  v_rel_after = Rotate(v_rel, v_rel_after1);

  vi = v_com + v_rel_after * b2;
  //velocity = v_com - v_rel_after * b1;8
  ion.velocity = vi;
  //std::cout << "after collision vi = " << vi << "\n";
}


void NeutralAtom::LoadThetaCache(unsigned int count) {
  double d = M_PI / count;
  theta_edge.resize(0);
  theta_center.resize(0);

  for (unsigned int i=0; i<count+1; i++) {
    theta_edge.push_back(i*d);
    theta_center.push_back((i+0.5) * d);
  }
  theta_center.pop_back();
}


void NeutralAtom::LoadVnCache(unsigned int count) {
  double d = 50.0 * sigma / count; 
  vn_edge.resize(0);
  vn_center.resize(0);

  for (unsigned int i=0; i<count+1; i++) {
    vn_edge.push_back(i*d);
    vn_center.push_back((i+0.5) * d);
  }
  vn_center.pop_back();
}

// void NeutralAtom::ScatteringTest(const string &fname) {
//   fstream f1((fname+".scattering").c_str(), fstream::out);
//   if (!f1) throw "File open wrong";
  
//   Vector3D vi_orig(-1, 3, 2);
//   Vector3D vi = vi_orig;
//   const int count = 100000;
//   double mi = 173*SI::AMU;
//   double qi = SI::ELECTRON_CHARGE;

//   for (int i=0; i<count; i++) {
//     Collide(mi, vi);
//   }
  
//   f1.close();
// }


// void NeutralAtom::SampleVelocityTest(const string &fname) {
//   fstream f1((fname+".vel").c_str(), fstream::out);
//   if (!f1) throw "File open wrong";

//   fstream f2((fname+".summary").c_str(), fstream::out);
//   if (!f2) throw "File open wrong";
  
  
//   Vector3D vi(-1, 3, 2);
//   const int count = 100000;
//   double rate_mc=0, rate=0;
//   double rate_integration=0;
//   double mi = 173*SI::AMU;
//   double qi = SI::ELECTRON_CHARGE;

//   for (int i=0; i<count; i++) {
//     SampleVelocity(vi);
//     f1 << velocity << "\n";
//     rate = GetCollisionRate(mi, qi, vi);
//     rate_mc += rate;
//   }
  
//   rate_mc /= count;
//   rate_integration = GetThermalAverageCollisionRate(mi, qi, vi);
//   f2 << "Rate Coeff from Monte Carlo = " << rate_mc << "\n";
//   f2 << "Rate Coeff from Integration = " << rate_integration << "\n";
//   f1.close();
//   f2.close();
// }

void NeutralAtom::LoadCrossSection(LuaConfigReader &reader) {
  typedef vector<string> string_list;

  string field("cross_section.lookup_table");
  string_list fname_list = reader.GetValue< string_list >(field.c_str());

  for (const string& fname:fname_list) {
    cout << "Loading diff x section " << fname << "\n";
    double energy;

    double t, t_prev, s, s_prev;
    vector<double> t_array, w_array;

    ifstream f(fname.c_str());
    if (!f) throw "File open wrong";
   
    f >> energy;
    while (f >> t >> s) {
      t_array.push_back(t);

      if (t_array.size()>1) 
    	w_array.push_back( (s+s_prev)/2 * sin((t+t_prev)/2) * (t-t_prev) );

      t_prev = t;
      s_prev = s;
    }

    diff_cross_section_lktb[energy] = pc_dist(begin(t_array), 
					      end(t_array), 
					      begin(w_array));
    
    f.close();
  }
}

NeutralAtom::pc_dist &NeutralAtom::SelectThetaDist(double energy) {
  auto it = diff_cross_section_lktb.upper_bound(energy);
  if (it==diff_cross_section_lktb.end()) it--;
  //  std::cout << "using E = " << it->first << "\n";
  return it->second;
}

void NeutralAtom::LoadFLktb() {
  for (const double *ff : f)
    f_lktb[ff[0]] = ff[1];
}

// double NeutralAtom::GetCollisionRate(double mi, double qi, const Vector3D &vi) {
//   double mu = mi * mass / (mi+mass);
//   double C4 = polarizability * qi * qi / (4*M_PI*EPSILON_0);
//   Vector3D v_rel = velocity - vi;

//   double energy = 0.5 * mu * v_rel.normSquared();
//   return v_rel.norm() * M_PI*pow(mu*pow(C4/HBAR,2)/energy, 1.0/3)*(1+pow(M_PI/4,2)) * density;
// }

