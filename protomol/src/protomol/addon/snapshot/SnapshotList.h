#ifndef _SNAPSHOT_LIST_H
#define _SNAPSHOT_LIST_H

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <protomol/addon/snapshot/TimeQueue.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/addon/snapshot/Snapshot.h>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <vector>
#include <string>

namespace ProtoMolAddon {
  namespace Snapshot {

    namespace pt = boost::property_tree;
    namespace algorithm = boost::algorithm;
    using namespace ProtoMol;
    
    template <class Storage> 
    class SnapshotList {

    private:
      std::vector< Snapshot<Storage> > ssl;
      static const std::string default_pattern;
      
    public:
      SnapshotList() {}
      
      SnapshotList(const std::string &fname) {
	pt::ptree tree;
 	pt::read_xml(fname, tree);

	std::string root = "ConfigRoot.SnapshotList" + Storage::GetName();
	std::string fname_pattern = tree.get<std::string>(root+".FileNamePattern",
							  default_pattern);
	algorithm::trim(fname_pattern);
	Snapshot<Storage>::SetFnamePattern(fname_pattern);
	
	for (auto &v: tree.get_child(root)) {
	  cout << v.first << "\n";
	  if (v.first=="Snapshot") ssl.push_back(Snapshot<Storage>(v.second));
	}
      }

      void Initialize(const ProtoMolApp *app) {
	for (auto &ss: ssl) ss.Initialize(app);
      }

      void Update(double now) {
	for (auto &ss: ssl) ss.Update(now);
      };

      void Finalize() {
	for (auto &ss: ssl) ss.Finalize();
      }

      static std::string GetName() { return "SnapshotList"+Storage::GetName(); }
      static std::string GetParameterName() { return "-snapshotlist"+Storage::GetParameterName(); }
      
    };

    template<class Storage>
    const std::string SnapshotList<Storage>::default_pattern("snapshot_%d.h5");
  }
}
//     private:
//       std::vector<Storage> storage_list;

//     public:
//       Snapshot() : storage_list(0) {}

//       Snapshot(const std::string &fname) {
// 	pt::ptree tree;
// 	pt::read_xml(fname, tree);

// 	std::string root = "ConfigRoot.OutputSnapshotWith" + Storage::GetName();
// 	boost::optional<std::string> pattern = tree.get_optional<std::string>(root+".FileNamePattern");
	
// 	if (pattern) 
// 	  Storage::SetFileNamePattern(pattern.get());

// 	for(auto &v: tree.get_child(root+".TimeQueue")) 
// 	  if (v.first=="Entry")
// 	    storage_list.push_back(Storage(v.second.get<TimeQueue>("")));
//       }
      

//       void Initialize(const ProtoMolApp *app) {
// 	for (auto &s: storage_list)
// 	  s.Initialize(app);
//       }

//       void Finalize() {
// 	for (auto &s: storage_list)
// 	  s.Finalize();
//       }

//       void Run(double now) {
// 	for (auto &s: storage_list)
// 	  s.SaveFrame(now);
//       }

//     };

//   }
// }
 
#endif 
