#ifndef _GENERIC_STORAGE_H
#define _GENERIC_STORAGE_H

#include <string>

using std::string;

namespace ProtoMol {
  class ProtoMolApp;
}

namespace ProtoMolAddon {
  namespace Snapshot {
    
    class GenericStorage {
      static unsigned int counter;
      static string fname_pattern;

      static void SetFilenamePattern(const string &pattern);
      
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
