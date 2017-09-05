#ifndef _TMDEVOL_H_
#define _TMDEVOL_H_

#include <iostream>
#include <fstream>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

#include "TMath.h"
#include "Math/GSLIntegrator.h"
#include "Math/SpecFuncMathMore.h"

using namespace std;

namespace TMDEVOL{

  const LHAPDF::PDF * xpdf = LHAPDF::mkPDF("CJ15lo", 0);

  const double gammaE = TMath::EulerGamma();
  const double Pi = TMath::Pi();

  const double C_F = 4.0 / 3.0;
  const double C_A = 3.0;
  const double T_R = 0.5;

  const double M_charm = 1.3;
  const double M_bottom = 4.5;
  const double M_top = 180.0;

  int order_C;
  int order_A;
  int order_B;
  
  double C1 = 2.0 * exp(-gammaE);

  double bstar_CS(const double b_T){//Collins and Soper, NPB 197 (1982) 446
    double bmax = 0.5;
    return b_T / sqrt(1.0 + pow(b_T / bmax, 2));
  }

  double bstar_bT(const double b_T){
    double bmax = 0.5;
    if (b_T < bmax) return b_T;
    else return bmax;
  }

  double (* bstar)(const double b_T);

  double mu_b(const double b_T){
    return C1 / bstar(b_T);
  }

  double N_f(const double mu){
    double value = 3.0;
    if (mu > M_charm) value = 4.0;
    if (mu > M_bottom) value = 5.0;
    if (mu > M_top) value = 6.0;
    return value;
  }
  
  double Afactor(const double mu){
    double value = 0;
    if (order_A > 0){
      value += (xpdf->alphasQ(mu) / Pi) * C_F;
    }
    if (order_A > 1){
      value += pow(xpdf->alphasQ(mu) / Pi, 2) * C_F / 2.0 * (C_A * (67.0 / 18.0 - Pi * Pi / 6.0) - 10.0 / 9.0 * T_R * N_f(mu));
    }
    return value;
  }

  double Bfactor(const double mu){
    double value = 0;
    if (order_B > 0){
      value += (xpdf->alphasQ(mu) / Pi) * (-3.0 * C_F / 2.0);
    }
    if (order_B > 1){
      double zeta3 = 1.202;//riemann_zeta(3.0)
      value += pow(xpdf->alphaQ(mu) / Pi, 2) * (pow(C_F / 2.0, 2) * (Pi * Pi - 3.0 / 4.0 - 12.0 * zeta3) + C_F / 2.0 * C_A * (11.0 / 18.0 * Pi * Pi - 193.0 / 24.0 + 3.0 * zeta3) + C_F / 2.0 * T_R * N_f(mu) * (17.0 / 6.0 - 2.0 / 9.0 * Pi * Pi));
    }
    return value;
  }

  double S_integrand(const double logmu, void * par){
    double Q = ((double *) par)[0];
    double mu = exp(logmu);
    double result = 2.0 * log(Q / mu) * A(mu) + B(mu);
    return result;
  }

  double S(const double b_T, const double Q){
    double par = Q;
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-6);
    ig.SetFunction(&S_integrand, &par);
    double result = ig.Integral(log(C1 / b_T), log(Q));
    return result;
  }

  double C_ij_1(const int flavor_i, const int flavor_j, const double mu){
    double value = 0.0;
    if (flavor_i == flavor_j) value = 1.0;
    if (order_C > 0){
      if (flavor_i == flavor_j){
	value += (xpdf->alphasQ(mu) / Pi) * C_F / 2.0 * (Pi * Pi / 2.0 - 4.0);
      }
    }
    return value;
  }

  double C_ij(const int flavor_i, const int flavor_j, const double mu, const double xi){
    double value = 0.0;
    if (order_C > 0){
      if (flavor_i == flavor_j){
	value += (xpdf->alphaQ(mu) / Pi) * C_F / 2.0 * (1.0 - xi);
      }
    }
    return value;
  }

  double C_ig_1(const int flavor_i, const double mu){
    double value = 0.0;
    if (order_C > 0.0){
      value += 0.0;
    }
    return value;
  }

  double C_ig(const int flavor_i, const double mu, const double xi){
    double value = 0.0;
    if (order_C > 0){
      value += (xpdf->alphasQ(mu) / Pi) * T_R * xi * (1.0 - xi);
    }
    return value;
  }

  double (* F_input)(const int flavor, const double x, const double b_T);


  double Q0_model0;
  double kt_model0;
  double F_input_model0(const int flavor, const double x, const double b_T){
    double Q0 = Q0_model0;
    double kt = model0;
    double S0 = S(b_T, Q0);
    double result = xpdf->xfxQ(flavor, x, Q0) / x * exp(-0.25 * pow(kt * b_T, 2));
    return result * exp(S0);
  }

  double F_output(const int flavor, const double x, const double b_T, const double Q){
    return F_input(flavor, x, b_T) * exp(-S(b_T, Q));
  }		      		    

  int Initialize(){
    order_C = 0;
    order_A = 2;
    order_B = 1;
    bstar = & bstar_bT;
    F_input = & F_input_model0;
    return 0;
  }
}





#endif
