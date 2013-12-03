#ifndef _ELECTRODE_H
#define _ELECTRODE_H
#define BOOST_DISABLE_ASSERTS
#include <protomol/type/Vector3D.h>
#include "boost/multi_array.hpp"
#include <vector>
#include <string>
#include <bitset>
#include <iostream>
#include <cassert>
#include <cmath>

using namespace std;
using namespace boost;
using namespace ProtoMol;


namespace ProtoMolAddon {

  class Electrode {
  public:
    class Voltage {
      vector<double>* v;
      double t0, dt;
      bool is_mirror;

    public:
      Voltage(): v(NULL), t0(0), dt(0) {}

      Voltage(vector<double>* v, double t0, double dt): v(v), t0(t0), dt(dt), is_mirror(false) {}

      Voltage& operator= (const Voltage& other) { 
	v = other.v;
	t0 = other.t0;
	dt = other.dt;
	is_mirror = true;
	return (*this);
      }

      ~Voltage() { if (v && !is_mirror) delete v; }

      bool Initialized() { return v!=NULL; }

      double GetVoltage(double t, int offset=0) {
	assert(v);
	int n = static_cast<int>((t - t0)/dt) + offset;
	
	if (n>v->size()-1)
	  n = v->size()-1;
	else if (n<0)
	  n = 0;
	
	//if (n>1000)
	//  cout << "t = " << t << "\tn = " << n << "\tv = " << v->at(n) << "\n";
	return v->at(n);
      }

      template <class Iterator>
      Voltage& Assign(double t0, double dt, int size, Iterator begin, Iterator end) {
	this->t0 = t0;
	this->dt = dt;
	if (v)
	  v->assign(begin, end);
	else 
	  v = new vector<double>(begin, end);

	//cout << v->size() << "\t" << size << "\n";
	assert(v->size() == size);
	return (*this);
      }

      friend istream& operator>> (ifstream& is, Voltage& volt) {
	double t0, dt;
	int size;
	is >> t0 >> dt >> size;
	istream_iterator<double> begin(is), end;
	volt.Assign(t0, dt, size, begin, end);
	return is;
      }


    };

    class Potential {
    private:
      multi_array<double, 3>* p;
      Vector3D x0, dx;
      bitset<3> reflection;
      bool is_mirror;

    public:
      Potential(): p(NULL), x0(Vector3D()), dx(Vector3D()), is_mirror(false) {};
      
      Potential(multi_array<double, 3>* p, const Vector3D& x0, const Vector3D& dx): p(p), x0(x0), dx(dx), reflection("") {}

      const Vector3D& GetDx() { return dx; }      
      Potential& operator= (const Potential& other) {
	p = other.p;
	x0 = other.x0;
	dx = other.dx;
	reflection = other.reflection;
	is_mirror = true;
	return (*this);
      }

      ~Potential() { if (p && !is_mirror) delete p; }

      bool Initialized() { return p!=NULL; }

      template <class Iterator>
      Potential& Assign(const Vector3D& x0, const Vector3D& dx, const boost::array<size_t, 3>& size, Iterator begin, Iterator end) {
	
	this->x0 = x0;
	this->dx = dx;
	if (p) 
	  p->resize(boost::extents[size[0]][size[1]][size[2]]);
	else 
	  p = new multi_array<double, 3>(boost::extents[size[0]][size[1]][size[2]]);
	
	p->assign(begin, end);	       
	assert(p->shape()[0] == size[0] && p->shape()[1] == size[1] && p->shape()[2] == size[2]);

	return (*this);
      }

      friend istream& operator>> (ifstream& is, Potential& potl) {
	Vector3D x0, dx;
	boost::array<size_t, 3> size;
	is >> x0 >> dx >> size[0] >> size[1] >> size[2];
	istream_iterator<double> begin(is), end;
	potl.Assign(x0, dx, size, begin, end);
	return is;
      }

      double GetPotential(const Vector3D& pos, const boost::array<int, 3>& offset) {
	assert(p!=NULL);
	typedef boost::multi_array<double, 3>::index index;
	boost::array<index, 3> n;
	for (int i=0; i<3; i++) {
	  int f = -2 * reflection[2-i] + 1;
	  n[i] = static_cast<int>((pos[i]*f-x0[i])/dx[i]) + offset[i]*f + reflection[2-i];
	  if (n[i] < 0 || n[i] > p->shape()[i]-1 ) return 0; 
	}
	return (*p)(n);
      }
      
      double GetInterpolatedPotential(const Vector3D& pos) {
	boost::array<double, 3> f;
	for (int i=0; i<3; i++) 
	  f[i] = fmod(pos[i]-x0[i], dx[i]);

	vector< bitset<3> > offset;
	for (int i=0; i<8; i++) 
	  offset.push_back(i);

	vector<double> weight;
	transform(offset.begin(), offset.end(), back_inserter(weight), [&f](const bitset<3>& o) {
	    double w = 1;
	    for (int i=0; i<3; i++) w *= (o[i] == 0 ? 1-f[i] : f[i]);
	    return w;
	  });

	vector<double> p;
	transform(offset.begin(), offset.end(), back_inserter(p), [&pos, this](const bitset<3>& o) { 
	    boost::array<int, 3> oi;
	    for (int i=0; i<3; i++) oi[i]=o[i];
	    return GetPotential(pos, oi); 
	  });

	return inner_product(p.begin(), p.end(), weight.begin(), 0.0);
      }

      void SetReflection(const string& r) {
	reflection = bitset<3>(r);
      };
    };

  private:
    string label;
    Voltage volt;
    Potential potl;
 
  public:
    Electrode(const string& label = "electrode") :
      label(label),
      volt(),
      potl()
    {}

    ~Electrode() {};
    void SetLabel(const string& label) {this->label = label; }
    const string& GetLabel() { return label; }

    Potential& Potl() { return potl; }
    Potential& Potl(const Potential& p) { potl = p; return potl;}
    
    Voltage& Volt() { return volt; }
    Voltage& Volt(const Voltage& v) { volt = v; return volt; }

    double GetRealTimePotential(const Vector3D& pos, double t, const boost::array<int, 3>& offset) {
      //      cout << label << "\n";
      return potl.GetPotential(pos, offset) * volt.GetVoltage(t);
    }

    double GetRealTimeInterpolatedPotential(const Vector3D& pos, double t) {
      return potl.GetInterpolatedPotential(pos) * volt.GetVoltage(t);
    }

    //template <class Iterator>
    //Electrode& Potl(const Vector3D& x0, const Vector3D& dx, array<int, 3>& ext, Iterator begin, Iterator end);

    //Electrode& MirrorPotl(const Electrode* p, const bitset<3> reflection);
    //Electrode& CopyVoltl(const Electrode* p);

    //double GetPotl(const Vector3D& pos, const vector<int, 3> offset);
    //double GetVolt(double t, int offset);
    friend ostream& operator<<(ostream& os, const Electrode& e) {
      os << "Electrode label = " << e.label << "\n";
      return os;
    }
  };

  
}

#endif
