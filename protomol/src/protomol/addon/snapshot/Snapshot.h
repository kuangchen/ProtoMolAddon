#ifndef _SNAPSHOT_H
#define _SNAPSHOT_H

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <protomol/addon/snapshot/TimeQueue.h>
#include <protomol/ProtoMolApp.h>
#include <fstream>
#include <vector>
#include <string>

using ProtoMol::ProtoMolApp;

namespace pt = boost::property_tree;

namespace ProtoMolAddon {
  namespace Snapshot {

    template <class Storage> 
    class Snapshot {
    private:
      std::vector<Storage> storage_list;

    public:
      Snapshot() : storage_list(0) {}

      Snapshot(const std::string &fname) {
	std::ifstream is(fname);
	if (!is)
	  throw std::runtime_error("Cannot open file " + fname);

	pt::ptree tree;
	pt::read_xml(fname, tree);

	std::string root = "SnapshotWith";
	
	for (auto &v: tree.get_child((root+Storage::GetName()).c_str())) {
	  if (v.first == Storage::GetName() + "Spec") {
	    Storage::SetFileNamePattern(v.second.get<std::string>("FileNamePattern"));
	    Storage::SetFrameNamePattern(v.second.get<std::string>("FileNamePattern"));
	  }

	  if (v.first == "TimeQueueSpec") 
	    storage_list.push_back(Storage(v.second));
	}
      }

      void Initialize(const ProtoMolApp *app) {
	for (auto &s: storage_list) s.Initialize(app);
      }

      void Finalize() {
	for (auto &s: storage_list) s.Finalize();
      }

      void Run(double now) {
	for (auto &s: storage_list) s.SaveFrame(now);
      }

    };

  }
}
 
#endif 
