#ifndef _SNAPSHOT_H
#define _SNAPSHOT_H

#include <boost/format.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <protomol/ProtoMolApp.h>
#include <string>
#include <memory>
#include <vector>
#include <queue>

namespace ProtoMol {
  class ProtoMolApp;
}

namespace ProtoMolAddon {

  namespace Util {
    class ConstSIAtomProxyArray;
  }
  
  namespace Snapshot {
    namespace pt = boost::property_tree;
    
    template <class Storage>
    class Snapshot {
    private:
      static size_t counter;
      static std::string pattern;
      
      Storage storage;
      std::queue<double> unsaved_time;
      std::vector<double> saved_time;
      
    public:
      Snapshot() {}

      Snapshot(const pt::ptree &tree) :
	storage((boost::format(pattern) % counter++).str())
      {
	double t0 = tree.get<double>("t0");
	size_t n = tree.get<size_t>("n");
	double dt = tree.get<double>("dt");

	for (unsigned int i=0; i<n; i++) unsaved_time.push(t0+i*dt);
      }

      Snapshot(const Snapshot &other):
	storage(other.storage),
	unsaved_time(other.unsaved_time),
	saved_time(other.saved_time)
      {}

      static void SetFnamePattern(const std::string &fname_pattern) {
	pattern = fname_pattern ;
      }
	
      void Initialize(const ProtoMol::ProtoMolApp *app) {
	storage.Initialize(app, unsaved_time);
      }

      void Update(double now) {
	// Don't do anything if everything has been saved
	// or the first element of unsaved queue is bigger than now
	if (unsaved_time.empty() || unsaved_time.front() > now) return;

	//
	storage.Save(unsaved_time, saved_time);
	//
	double t = unsaved_time.front();
	unsaved_time.pop();
	saved_time.push_back(t);
      }
      
      void Finalize() {
	storage.Finalize(unsaved_time, saved_time);
      }
    };

    template<class Storage>
    size_t Snapshot<Storage>::counter = 0;

    template<class Storage>
    std::string Snapshot<Storage>::pattern("snapshot_%d.h5");
  }
}


//     namespace pt = boost::property_tree;
//     using namespace ProtoMol;
    
//     template <class Storage> 
//     class Snapshot {
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
