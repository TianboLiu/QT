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

  LHAPDF::PDF * xpdf;

  const double gammaE = TMath::EulerGamma();//Euler 
  const double Pi = TMath::Pi();//pi

  const double C_F = 4.0 / 3.0;//color factors
  const double C_A = 3.0;
  const double T_R = 0.5;

  int order_A = 2;//orders
  int order_B = 1;
  int order_C = 0;
  
  double C1 = 2.0 * exp(-gammaE);

  double (* bstar)(const double b_T);

  double bstar_CS(const double b_T){//Collins and Soper, NPB 197 (1982) 446
    double bmax = 0.5;
    return b_T / sqrt(1.0 + pow(b_T / bmax, 2));
  }

  double bstar_sharp(const double b_T){
    double bmax = 0.5;
    if (b_T < bmax) return b_T;
    else return bmax;
  }

  double bstar_bT(const double b_T){
    return b_T;
  }
   
  double mu_b(const double b_T){
    return C1 / bstar(b_T);
  }

  double N_f(const double mu){
    double value = 3.0;
    if (mu > xpdf->quarkThreshold(4)) value = 4.0;
    if (mu > xpdf->quarkThreshold(5)) value = 5.0;
    if (mu > xpdf->quarkThreshold(6)) value = 6.0;
    return value;
  }
  
  double Afactor(const double mu){
    double value = 0;
    // order 0
    if (order_A < 1) return value;
    // order 1
    value += (xpdf->alphasQ(mu) / Pi) * C_F;
    if (order_A < 2) return value;
    // order 2
    value += pow(xpdf->alphasQ(mu) / Pi, 2) * C_F / 2.0 * (C_A * (67.0 / 18.0 - Pi * Pi / 6.0) - 10.0 / 9.0 * T_R * N_f(mu));
    if (order_A < 3) return value;
    return value;
  }

  double Bfactor(const double mu){
    double value = 0;
    //order 0
    if (order_B < 1) return value;
    //order 1
    value += (xpdf->alphasQ(mu) / Pi) * (-3.0 * C_F / 2.0);
    if (order_B < 2) return value;
    //order 2
    double zeta3 = 1.202;//riemann_zeta(3.0)
    value += pow(xpdf->alphasQ(mu) / Pi, 2) * (pow(C_F / 2.0, 2) * (Pi * Pi - 3.0 / 4.0 - 12.0 * zeta3) + C_F / 2.0 * C_A * (11.0 / 18.0 * Pi * Pi - 193.0 / 24.0 + 3.0 * zeta3) + C_F / 2.0 * T_R * N_f(mu) * (17.0 / 6.0 - 2.0 / 9.0 * Pi * Pi));
    if (order_B < 3) return value;
    return value;
  }

  double S_integrand(const double logmu, void * par){
    double Q = ((double *) par)[0];
    double mu = exp(logmu);
    double result = 2.0 * log(Q / mu) * Afactor(mu) + Bfactor(mu);
    return result;
  }

  double S(const double b_T, const double Q){
    double par = Q;
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-6);
    ig.SetFunction(&S_integrand, &par);
    double result = ig.Integral(log(mu_b(b_T)), log(Q));
    return result;
  }

  double C_ij_1(const int i, const int j, const double mu){
    //delta(z - 1) part
    double value = 0.0;
    //order 0
    if (i == j){
      value += 1.0;
    }
    if (order_C < 1) return value;
    //order 1
    if (i == j && i != 21){
      value += (xpdf->alphasQ(mu) / Pi) * C_F / 2.0 * (Pi * Pi / 2.0 - 4.0);
    }
    if (order_C < 2) return value;
    return value;
  }

  double C_ij(const int i, const int j, const double z, const double mu){
    //non-delta part
    double value = 0.0;
    //order 0
    if (order_C < 1) return value;
    //order 1
    if (i == j && j != 21){
      value += (xpdf->alphasQ(mu) / Pi) * C_F / 2.0 * (1.0 - z);
    }
    if (i != 21 && j == 21){
      value += (xpdf->alphasQ(mu) / Pi) * T_R * z * (1.0 - z);
    }
    if (order_C < 2) return value;
    return value;
  }

  double F_col_integrand(const double xi, void * par){
    //par: flavor, x, mu
    double * para = (double *) par;
    int flavor = floor(para[0]);
    double x = para[1];
    double mu = para[2];
    double z = x / xi;
    double sum = 0.0;
    for (int j = 1; j <= 6; j++){
      sum += C_ij(flavor, j, z, mu) * xpdf->xfxQ(j, xi, mu);
    }
    for (int j = -1; j >= -6; j--){
      sum += C_ij(flavor, j, z, mu) * xpdf->xfxQ(j, xi, mu);
    }
    sum += C_ij(flavor, 21, z, mu) * xpdf->xfxQ(21, xi, mu);
    return sum / xi;
  }

  double F_col(const int flavor, const double x, const double b_T){
    double par[3] = {((double) flavor) + 1e-3, x, mu_b(b_T)};
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-6);
    ig.SetFunction(&F_col_integrand, par);
    double value = ig.Integral(x, 1.0 - 1e-9);
    for (int j = 1; j <= 6; j++){
      value += C_ij_1(flavor, j, mu_b(b_T)) * xpdf->xfxQ(j, x, mu_b(b_T));
    }
    for (int j = -1; j >= -6; j--){
      value += C_ij_1(flavor, j, mu_b(b_T)) * xpdf->xfxQ(j, x, mu_b(b_T));
    }
    value += C_ij_1(flavor, 21, mu_b(b_T)) * xpdf->xfxQ(21, x, mu_b(b_T));
    return value;
  }    
  
  double (* g_K)(const double b_T);

  double g_K_zero(const double b_T){
    return 0.0;
  }
  
  double (* F_NP)(const int flavor, const double x, const double b_T);

  double F_NP_gaus(const int flavor, const double x, const double b_T){
    return exp(-0.25 * 0.5 * b_T * b_T);
  }

  double Q0 = 1.0;//input scale

  double F_TMD(const int flavor, const double x, const double b_T, const double Q){
    double value = F_col(flavor, x, b_T) * exp(-S(b_T, Q)) * exp(g_K(b_T) * log(Q / Q0)) * F_NP(flavor, x, b_T);
    return value;
  }		      		    

  int Initialize(){
    order_C = 0;
    order_A = 2;
    order_B = 1;
    bstar = & bstar_bT;
    g_K = & g_K_zero;
    F_NP = & F_NP_gaus;
    xpdf = LHAPDF::mkPDF("CJ15lo", 0);
    return 0;
  }
}





#endif
