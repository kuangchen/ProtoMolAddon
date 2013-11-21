#ifndef _ELECTRODE_H
#define _ELECTRODE_H

#include <protomol/type/Vector3D.h>
#include "boost/multi_array.hpp"
#include <vector>
#include <string>
#include <bitset>
#include <iostream>
#include <cassert>


using namespace std;
using namespace boost;
using namespace ProtoMol;


namespace ProtoMolAddon {

  class Electrode {
  public:
    class Voltage {
      vector<double>* v;
      double t0, dt;

    public:
      Voltage(): v(NULL), t0(0), dt(0) {}

      Voltage(vector<double>* v, double t0, double dt): v(v), t0(t0), dt(dt) {}

      Voltage& operator= (const Voltage& other) { 
	v = other.v;
	t0 = other.t0;
	dt = other.dt;
	return (*this);
      }

      ~Voltage() { if (v) delete v; }

      double GetVolt(double t, int offset=0) {
	assert(v);
	int n = static_cast<int>((t - t0)/dt) + offset;
	
	if (n>v->size())
	  n = v->size();
	else if (n<0)
	  n = 0;

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
	
	assert(v->size() == size);
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

    public:
      Potential() {};
      
      Potential(multi_array<double, 3>* p, const Vector3D& x0, const Vector3D& dx): p(p), x0(x0), dx(dx), reflection("") {}
      
      Potential& operator= (const Potential& other) {
	p = other.p;
	x0 = other.x0;
	dx = other.dx;
	reflection = other.reflection;
	return (*this);
      }

      ~Potential() { if (p) delete p; }

      template <class Iterator>
      Potential& Assign(const Vector3D& x0, const Vector3D& dx, const boost::array<size_t, 3>& size, Iterator begin, Iterator end) {
	this->x0 = x0;
	this->dx = dx;
	if (p) 
	  p->resize(size);
	else 
	  p = new multi_array<double, 3>(size);
	
	p->assign(begin, end);	       
	assert(p->shape()[0] == size[0] && p->shape()[1] == size[1] && p->shape()[2] == size[2]);
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
	assert(p);
	boost::array<int, 3> n;
	for (int i=0; i<3; i++) {
	  n[i] = (static_cast<int>((pos[i]-x0[0])/dx[i]) + offset[i]) * reflection[2-i];
	  if (n[i] < 0 || n[i] > p->shape()[i]) 
	    return 0;
	}
	return (*p)(n);
      }
      
      void SetReflection(const bitset<3>& r) {
	reflection = r;
      };
    };

  private:
    string label;
    Voltage volt;
    Potential potl;

  public:
    Electrode(const string& label) :
      label(label),
      volt(),
      potl()
    {}

    ~Electrode() {};

    Potential& Potl() { return potl; }
    Potential& Potl(const Potential& Potl) { this->potl = potl; }
    
    Voltage& Volt() { return volt; }
    Voltage& Volt(const Voltage& Potl) { this->volt = volt; }
    
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
