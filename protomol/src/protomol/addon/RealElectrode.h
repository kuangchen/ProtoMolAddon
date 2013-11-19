#ifndef __REAL_ELECTRODE_H
#define __REAL_ELECTRODE_H

#include <protomol/addon/AbstractElectrode.h>
#include "boost/multi_array.hpp"
#include <string>
#include <vector>
#include <array>

using namespace ProtoMol;
using namespace std;

namespace ProtoMolAddon {
  class RealElectrode : public AbstractElectrode {
  private:
    vector<double> volt_data;
    double t0;
    double dt;

    boost::multi_array<double, 3> pot_data;
    Vector3D x0;
    Vector3D dx;
    
  public:
    typedef boost::multi_array<double, 3>::size_type size_type;
    typedef boost::array<size_type, 3> extent_type;
    typedef boost::multi_array<double, 3>::extent_gen extent_gen;
    typedef boost::multi_array<double, 3>::index index;
    typedef boost::array<index, 3> indices_type;

  public:
    RealElectrode(const string& label) : AbstractElectrode(label), volt_data(0), t0(0), dt(0), x0(), dx() {}
    ~RealElectrode() {}

    template<class Iterator>
    RealElectrode& AssignVoltage(double t0, double dt, Iterator begin, Iterator end) {
      this->t0 = t0;
      this->dt = dt;
      volt_data.assign(begin, end);
      return (*this);
    }

    RealElectrode& AssignVoltageFromFile(const string& volt_file) {
      ifstream file(volt_file);
      double t0_ = 0, dt_ = 0;
      file >> t0_ >> dt_;
      istream_iterator<double> begin(file), end;
      AssignVoltage(t0_, dt_, begin, end);
      file.close();
      return (*this);
    }

    RealElectrode& AssignPotentialFromFile(const string& pot_file) {
      ifstream file(pot_file);
      Vector3D x0_, dx_;
      boost::array<int, 3> ext;
      file >> x0_[0] >> x0_[1] >> x0_[2];
      file >> dx_[0] >> dx_[1] >> dx_[2];
      file >> ext[0] >> ext[1] >> ext[2];

      // cout << "label = " << label << "\n";
      // cout << x0_[0] << " " << x0[1] << " " << x0[2];
      // cout << dx_[0] << " " << dx_[1] << " " << dx_[2];
      // cout << ext[0] << " " << ext[1] << " " << ext[2];

      istream_iterator<double> begin(file), end;
      AssignPotential(x0_, dx_, ext, begin, end);
      file.close();
      return (*this);
    }

    template<class Iterator>
    RealElectrode& AssignPotential(const Vector3D& x0, 
				   const Vector3D& dx, 
				   boost::array<int, 3> ext, 
				   Iterator begin, 
				   Iterator end) {
      this->x0 = x0;
      this->dx = dx;
      
      pot_data.resize(boost::extents[ext[0]][ext[1]][ext[2]]);
      pot_data.assign(begin, end);
      return (*this);
    }

    void GetFraction(const Vector3D &pos, array<double, 3>& f) {
      array<int, 3> nn;
      for(int i=0; i<3; i++) {
	nn[i] = static_cast<int>( (pos[i]-x0[i]) / dx[i] );
	f[i] = (pos[i] - x0[i] - nn[i] * dx[i])/dx[i];
      }
    }

    double GetNNPotential(const Vector3D& pos, const array<int, 3>& offset) const {
      indices_type nn;
      for (int i=0; i<3; i++) 
	nn[i] = static_cast<int>( (pos[i]-x0[i]) / dx[i] ) + offset[i];
      
      return pot_data(nn);
    }

    double GetNNVoltage(double time, int offset=0) const {
      int nn = static_cast<int>((time - t0)/dt);

      if (nn+offset > volt_data.size())
	return volt_data[volt_data.size()];

      else return volt_data[nn+offset];
    }

    const Vector3D& GetDx() {
      return this->dx;
    }

    void DumpInfo(ostream& os) {
      os << label << "\n"
	 << "x0 = (" << x0 << ")\n" 
	 << "dx = (" << dx << ")\n" 
	 << "size = (" << pot_data.shape()[0] << "," << pot_data.shape()[1] << "," << pot_data.shape()[2] << ")\n"
	 << "t0 = " << t0
	 << "dt = " << dt
	 << "size = " << volt_data.size() << "\n";
    }
  };

}


#endif
