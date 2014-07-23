#include <protomol/addon/ionization/IonizationManager.h>
#include <cmath>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>

using namespace ProtoMolAddon::Constant;
using namespace ProtoMolAddon::Ionization;
namespace pt = boost::property_tree;
namespace algorithm = boost::algorithm;

IonizationManager::Spec::Spec(const std::string &fname) :
  entry_list() {
  
  pt::ptree tree;
  pt::read_xml(fname, tree);

  for (auto &v: tree.get_child("ConfigRoot.IonizationSpec"))
    if (v.first=="entry") {
      std::string name(v.second.get<std::string>("atom_name"));
      algorithm::trim(name);
      double rate = v.second.get<double>("ionization_rate");
      int charge = v.second.get<int>("charge");

      entry_list.push_back(entry(name, rate, charge));
    }
}

void IonizationManager::Init(ProtoMol::ProtoMolApp *app) {
  for (unsigned int i=0; i<app->positions.size(); i++) {
    atom_proxy_list.push_back(Util::SIAtomProxy(app, i));
    for (unsigned int j=0; j<spec.entry_list.size(); j++) {
      if (app->topology->atoms[i].name == spec.entry_list[j].atom_name) {
	proxy_map[i] = j;
	continue;
      }
    }
  }
}

void IonizationManager::Ionize(double dt) {
  for (auto &v: proxy_map) {
    Util::SIAtomProxy &atom_proxy = atom_proxy_list[v.first];
    Spec::entry &e = spec.entry_list[v.second];

    if (atom_proxy.GetCharge() == 0) {
      double rate = e.ionization_rate;
      int charge = e.charge;
      
      if (exp(-rate * dt) < dist(engine)) {
	std::cout << "Atom " << atom_proxy.GetID() << " is ionized to charge " << charge << " state" << std::endl;
	atom_proxy.SetCharge(charge * ToSI::charge);
      }
    }
  }
}


