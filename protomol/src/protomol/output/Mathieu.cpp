#include <cmath>
#include <iostream>
#include "Mathieu.h"

using namespace Mathieu;
using std::cout;

MathieuFunction::MathieuFunction(double q, double a, int n_max):
  q(q), 
  a(a),
  n_max(n_max) 
{
  ksi.resize(n_max);
  alpha.resize(n_max);
  d.resize(n_max);

  coeff_n.resize(n_max+1);
  coeff_p.resize(n_max+1);

  ksi_mu_n.resize(n_max);
  ksi_mu_p.resize(n_max);

  ratio_n.resize(n_max);
  ratio_p.resize(n_max);
  
  int i=0;
  for (i = 0; i < n_max ; i++ ) 
    ksi[i] = q / (4 * i * i - a);

  // alpha[i] = ksi[i] * ksi[i-1]
  for (i = 1; i < n_max ; i++ )
    alpha[i] = ksi[i] * ksi[i-1];

  d[0] = 1;
  d[1] = 1 - 2 * alpha[1];
  d[2] = (2 * alpha[1] + alpha[2] - 1) * (alpha[2] - 1 );

  for (i = 3 ; i < n_max; i++ )
    d[i] = (1 - alpha[i]) * d[i-1] - alpha[i] * (1-alpha[i]) * d[i-2] + alpha[i] * alpha[i-1] * alpha[i-1] * d[i-3];

  double (*f)(double) = a >= 0 ? &cos : &cosh;
  double aa = fabs(a);

  mu = acos( 1 - d[n_max-1] * ( 1 - f(M_PI* sqrt(aa)))) / M_PI;

  for (i = 0; i < n_max ; i++ ) {
    double arg_p = (-2*i - mu);
    double arg_n = ( 2*i - mu);
    ksi_mu_n[i] = ( arg_p * arg_p - a ) / q;
    ksi_mu_p[i] = ( arg_n * arg_n - a ) / q;
  }

  ratio_n[n_max-1] = 0;
  ratio_p[n_max-1] = 0;

  for ( i = n_max-1 ; i > 0 ; i-- ) {
    ratio_n[i-1] = - 1 / (ratio_n[i] + ksi_mu_n[i]);
    ratio_p[i-1] = - 1 / (ratio_p[i] + ksi_mu_p[i]);
  }

  coeff_p[0] = 1;
  coeff_n[0] = 1;

  nf = 1;
  for ( i = 0; i < n_max-1 ; i++ ) {
    coeff_p[i+1] = coeff_p[i] * ratio_p[i];
    coeff_n[i+1] = coeff_n[i] * ratio_n[i];
    nf += (coeff_n[i+1] * coeff_n[i+1] + coeff_p[i+1] * coeff_p[i+1]);
  }

  nf = sqrt(nf);

  for ( i = 0; i < n_max ; i++ ) {
    coeff_p[i] /= nf;
    coeff_n[i] /= nf;
  }

  CalculateMiscVar (1e-8);
}

double MathieuFunction::GetMu() const {
  return mu;
}

void MathieuFunction::Evaluate(double x, vector<double>& r) {
  r.resize(0);
  int i;
  double arg = mu*x;
  double c = cos(arg) * coeff_p[0];
  double cp = -sin(arg) * mu * coeff_p[0];
  double s = sin(arg) * coeff_n[0];
  double sp = cos(arg) * mu * coeff_n[0];
  for ( i=1 ; i < n_max; i++ ) {
    double arg_p = mu - i*2.0;
    double arg_n = mu + i*2.0;

    double cos_p = cos(arg_p*x) * coeff_p[i];
    double cos_n = cos(arg_n*x) * coeff_n[i]; 

    double sin_p = sin(arg_p*x) * coeff_p[i];
    double sin_n = sin(arg_n*x) * coeff_n[i]; 

    c += cos_p;
    c += cos_n;

    s += sin_p;
    s += sin_n;

    cp -= sin_p * arg_p;
    cp -= sin_n * arg_n;

    sp += cos_p * arg_p;
    sp += cos_n * arg_n;
  }

  r[0] = c;
  r[1] = s;
  r[2] = cp;
  r[3] = sp;


}

void MathieuFunction::GetMiscVar(MathieuMiscVar *m) const {
  *m = MathieuMiscVar(misc_var);
}


void MathieuFunction::CalculateMiscVar(double accuracy) {
  double coeff_pi, coeff_ni;
  vector<double> r;
  r.resize(4);
  Evaluate(0, r);                     // c, s, cp, sp   
  misc_var.wronskian = r[0]*r[3] - r[1]*r[2]; // c*sp - s*cp
  
  misc_var.c_s_sq_sum_ave = 1;                // from normalization
  misc_var.epsilon = 0;

  double zeroth_ave = coeff_p[0] * coeff_p[0] * mu * mu;
  double curr = zeroth_ave;
  double prev = 0;

  int i = 1;
  do {
    prev = curr;
    coeff_pi = coeff_p[i];
    coeff_ni = coeff_n[i];

    curr += coeff_pi * coeff_pi * (mu-2*i) * (mu-2*i);
    curr += coeff_ni * coeff_ni * (mu+2*i) * (mu+2*i);
    
    i++;
  } while (fabs(prev/curr-1) > accuracy);
  misc_var.cp_sp_sq_sum_ave = curr;
  misc_var.eta = zeroth_ave / curr;

  curr = 0;
  prev = 0;

  double tau;
  double temp;

  i = 1;
  double dt = 1e-1;
  do {
    prev = curr;
    tau = dt * i;
    Evaluate (tau, r);
    temp = r[0] * r[2] + r[1] * r[3];
    curr = (curr*i + temp*temp)/(i+1);
    i++;
  } while (fabs(prev/curr-1) > accuracy);
  misc_var.epsilon = curr / (misc_var.wronskian * misc_var.wronskian);
}
