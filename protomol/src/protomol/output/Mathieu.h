#ifndef __MATHIEU_H
#define __MATHIEU_H

#include <vector>

using std::vector;

namespace Mathieu {

  struct MathieuMiscVar {
    double wronskian;           // Wronskian
    double c_s_sq_sum_ave;      // c^2 + s^2 averaged over time
    double cp_sp_sq_sum_ave;    // cp^2 + sp^2 averaged over time
    double epsilon;             // (c*cp + s*sp)^2 / wronskian^2
    double eta;                 // fraction of secular energy to total energy
  };
  
  class MathieuFunction {
  private:
    double q, a;
    int n_max;
    double mu;
    double nf;
    vector<double> ksi;
    vector<double> d;
    vector<double> alpha;

    vector<double> ratio_n;
    vector<double> ratio_p;
    vector<double> coeff_p;
    vector<double> coeff_n;
    vector<double> ksi_mu_p;
    vector<double> ksi_mu_n;

    struct MathieuMiscVar misc_var;

  private:
    void CalculateMiscVar(double accuracy);

  public:
    MathieuFunction(double q, double a, int n_max);

    void Evaluate(double x, vector<double>& r);
    double GetMu() const;
    void GetMiscVar(MathieuMiscVar *m) const;
  };
}

#endif
