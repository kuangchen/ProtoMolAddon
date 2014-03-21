#ifndef _GENERIC_STORAGE_H
#define _GENERIC_STORAGE_H

#include <string>

namespace ProtoMol {
  class ProtoMol;
}

using std::string;

namespace ProtoMolAddon {
  namespace Snapshot {
    
    class GenericStorage {
    private:
      static unsigned int counter;
      static string fname_pattern;
      static void SetFilenamePattern(const string &pattern) { fname_pattern = pattern; }
      
    public:
      GenericStorage();
	
      virtual ~GenericStorage() = 0; 
      virtual void Save(const ProtoMol::ProtoMolApp *app) = 0;

      unsigned int id;
      string fname;
    };

  }
}

#endif
