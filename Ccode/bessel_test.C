#include <iostream>
#include <fstream>
#include <cmath>

#include "Math/GSLIntegrator.h"
#include "Math/SpecFuncMathMore.h"

using namespace std;


double Accurate(const double width){
  return 0.5 * width * exp(-width / 4.0);
}

double integrand(const double ATan_b, void * par){
  double * params = (double *) par;
  double width = params[0];
  double b = tan(ATan_b);
  double result = b * ROOT::Math::cyl_bessel_j(0, b) * exp(- b * b / width);
  return result / pow(cos(ATan_b), 2);
}

double Numerical(const double width){
  double w = width;
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1e-8, 10000000);
  ig.SetFunction(&integrand, &w);
  double result = ig.Integral(0.0, M_PI_2 - 0.1) + ig.Integral(M_PI_2 - 0.1, M_PI_2);
  return result;
}

int main(const int argc, const char * argv[]){

  double width;
  double ratio;
  for (int i = 0; i < 1000; i++){
    width = 0.5 * (i + 1);
    ratio = Numerical(width) / Accurate(width);
    printf("%.2E \t %.6E \n", width, ratio);
  }

  return 0;
}
