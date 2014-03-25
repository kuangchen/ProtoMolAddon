#ifndef _DAMPING_SPEC_H
#define _DAMPING_SPEC_H

#include <protomol/ProtoMolApp.h>
#include <string>
#include <iosfwd>
#include <vector>
#include <map>

using std::ostream;
using std::string;
using std::istream;
using std::map;
using namespace ProtoMol;

namespace ProtoMolAddon {
  namespace Damping {

    class DampingSpec {
    public:
      struct DampingSpecEntry {
	double t_start;
	double t_end;
	double alpha;

	DampingSpecEntry(double t_start=0, double t_end=0, double alpha=0);
	void CheckValidity();
	friend istream& operator>> (istream &is, DampingSpecEntry &e);
	friend ostream& operator<< (ostream &os, const DampingSpecEntry &e);
      };

    private:
      map<string, DampingSpecEntry> entry_map;
   
    public:
      DampingSpec();
      void Damp(ProtoMolApp *app, double now, double h) const;

      friend istream& operator>> (istream &is, DampingSpec &spec);
      friend ostream& operator<< (ostream &os, const DampingSpec &spec);
    };

  }
}

#endif
