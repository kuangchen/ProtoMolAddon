#ifndef BUILD_TOPOLOGY_FROM_TPR_H
#define BUILD_TOPOLOGY_FROM_TPR_H

#include <protomol/topology/ExclusionType.h>
#include <protomol/type/Vector3DBlock.h>

#include <vector>
#include <string>


namespace ProtoMol {
  class GenericTopology;

  //build topo from TPR file
  void buildTopologyFromTpr(GenericTopology *topo, 
                                Vector3DBlock &pos, Vector3DBlock &vel,
                                const string &fname);
  //struct for saving functions
  struct function{
    string name;
    vector<double> parameters;
  };

  //types
  enum ftype{BOND, ANGLE, PROPERDIH, RYCKAERTBELL, LJ14, CONSTRAINT, FTSIZE };

  //struct for saving force field indexes
  struct ffdata{
      string name;
      string longname;
      ftype type;
      unsigned int function;
      vector<unsigned int> atoms;
  };

  //pair of unsigned ints for saving exclusions
  struct exclusion_pair{
      exclusion_pair(const unsigned int a, const unsigned int b){
          atom1 = a;
          atom2 = b;
      }
      unsigned int atom1;
      unsigned int atom2;
  };

  //parse functions
  bool parse_iparams(function &func, void *ft, void *ip,
                        ostringstream &os);

  //atomic radius from lookup, choose from orig. and Bowman sets
  double atom_radius( std::string atom_type, int set );

}

#endif // BUILD_TOPOLOGY_FROM_TPR_H
