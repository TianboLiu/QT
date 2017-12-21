/* Bessel function */
#ifndef _BESSEL_H_
#define _BESSEL_H_

#include <cmath>
#include "Math/SpecFuncMathMore.h"
#include "gsl/gsl_sf_bessel.h"
#include "Math/GSLIntegrator.h"

double Factorial(const int n){
  if (n == 0 || n == 1) return 1.0;
  double result = (double) n;
  for (int i = n - 1; i > 1; i--){
    result *= (double) i;
  }
  return result;
} 

double BesselJ_GSL(int n, double x){
  return gsl_sf_bessel_Jn(n, x);
}

double BesselJ_ROOT(int n, double x){
  return ROOT::Math::cyl_bessel_j(n, x);
}

double integrand(const double t, void * par){
  double * para = (double *) par;
  return cos(para[0] * t - para[1] * sin(t));
}

double BesselJ_Integral(int n, double x){
  double par[2] = {(double) n, x};
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1e-3, 10000000);
  ig.SetFunction(&integrand, par);
  return (ig.Integral(0.0, M_PI - 0.1) + ig.Integral(M_PI - 0.1, M_PI)) / M_PI;
}

double BesselJ_Series(int n, double x){
  double sum = 0;
  if (n == 0) sum = 1;
  else sum = pow(x / 2.0, n) / Factorial(n);
  int m = 1;
  double err = sum;
  while (abs(err / sum) > 1e-6){
    err = pow(x / 2.0, 2 * m + n) / Factorial(m) / Factorial(m + n);
    if (m%2 == 0)
      sum = sum + err;
    else
      sum = sum - err;
    m++;
    if (m > 20) break;
  }
  return sum;
} 

double BesselJ(int n, double x){
  return 0;
}







#endif
