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
  /* Following the formalism in 
     S.M. Aybat and T.C. Rogers, PRD 83 (2011) 114042
   */

  const LHAPDF::PDF * xpdf = LHAPDF::mkPDF("CJ15lo", 0);

  const double gammaE = TMath::EulerGamma();
  const double Pi = TMath::Pi();

  const double C_F = 4.0 / 3.0;
  const double C_A = 3.0;
  const double T_R = 0.5;

  double M_charm = 1.3;
  double M_bottom = 4.5;
  double M_top = 180.0;

  int order_C;
  int order_gamma_K;
  int order_gamma_F;
  
  double C1 = 2.0 * exp(-gammaE);

  double bstar_CS(const double b_T){//Collins and Soper, NPB 197 (1982) 446
    double bmax = 0.5;
    return b_T / sqrt(1.0 + pow(b_T / bmax, 2));
  }

  double (* bstar)(const double b_T);

  double mu_b(const double b_T){
    return C1 / bstar(b_T);
  }

  double kernel_K(const double b_T, const double mu){
    double value = - (xpdf->alphasQ(mu)) * C_F / Pi * 2.0 * (log(b_T * mu) - log(2.0) + gammaE);
    return value;
  }

  double gamma_K(const double mu){
    double value = 0.0;
    if (order_gamma_K > 0){
      value += 2.0 * (xpdf->alphasQ(mu) / Pi) * C_F;
    }
    if (order_gamma_K > 1){
      double Nf = 3.0;
      if (mu > M_charm) Nf = 4.0;
      if (mu > M_bottom) Nf = 5.0;
      if (mu > M_top) Nf = 6.0;
      value += 2.0 * pow(xpdf->alphasQ(mu) / Pi, 2) * C_F / 2.0 * (C_A * (67.0 / 18.0 - Pi * Pi / 6.0) - 10.0 / 9.0 * T_R * Nf);
    }
    return value;
  }

  double gamma_F(const double mu, const double zeta_F){
    double value = 0.0;
    if (order_gamma_F > 0){
      value += (xpdf->alphasQ(mu) / Pi) * C_F * (3.0 / 2.0 - log(zeta_F / (mu * mu)));
    }
    if (order_gamma_F > 1){
      double Nf = 3.0;
      if (mu > M_charm) Nf = 4.0;
      if (mu > M_bottom) Nf = 5.0;
      if (mu > M_top) Nf = 6.0;
      value += pow(C_F / 2.0, 2) * (12.0 * ROOT::Math::riemann_zeta(3.0) + 3.0 / 4.0 - Pi * Pi) - C_F * C_A / 2.0 * (Pi * Pi * 11.0 / 18.0 - 193.0 / 24.0 + 3.0 * ROOT::Math::riemann_zeta(3.0)) - C_F * T_R * Nf / 2.0 * (17.0 / 6.0 - Pi * Pi * 2.0 / 9.0);
    }
    return value;
  }

  double logB_integrand_2(const double logmu, void * par){
    double zeta_F = ((double *) par)[0];
    double mu = exp(logmu);
    return gamma_F(mu, mu * mu) - log(sqrt(zeta_F) / mu) * gamma_K(mu);
  }

  double logB_2(const double b_T, const double mu, const double zeta_F){
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-6);
    double par[1] = {zeta_F};
    ig.SetFunction(&logB_integrand_2, par);
    double min_logmu = log(mu_b(b_T));
    double max_logmu = log(mu);
    double result = ig.Integral(min_logmu, max_logmu);
    return result;
  }

  double B_factor(const double b_T, const double mu, const double zeta_F){
    double logB1 = log(sqrt(zeta_F) / mu_b(b_T)) * kernel_K(bstar(b_T), mu_b(b_T));
    double logB2 = logB_2(b_T, mu, zeta_F);
    return exp(logB1 + logB2);
  }

  double A_factor(const int flavor, const double x, const double b_T){
    double value = xpdf->xfxQ(flavor, x, mu_b(b_T)) / x;
    if (order_C > 0){
      value += 0.0;
    }
    return value;
  }

  double (* C_factor)(const int flavor, const double x, const double b_T, const double zeta_F, const double zeta_F_0);

  double C_factor_AR(const int flavor, const double x, const double b_T, const double zeta_F, const double zeta_F_0){
    double x0 = 0.02;
    double g1 = 0.21;
    double g2 = 0.68;
    double g3 = -0.6;
    double logC = (g2 / 2.0 * log(sqrt(zeta_F / zeta_F_0)) + g1 * (0.5 + g3 * log(10.0 * x * x0 / (x0 + x)))) * pow(b_T, 2);
    return exp(-logC);
  }

  double C_factor_Model0(const int flavor, const double x, const double b_T, const double zeta_F, const double zeta_F_0){
    double Ar = (xpdf->xfxQ(flavor, x, sqrt(zeta_F))) / (xpdf->xfxQ(flavor, x, mu_b(b_T)));
    double g2 = 0.68;
    double kt2 = 0.5;
    double C = exp(-(g2 / 2.0 * log(sqrt(zeta_F / zeta_F_0)) + 0.25 * kt2) * pow(b_T, 2));
    return C * Ar;
  }
 

  double F_bspace(const int flavor, const double x, const double b_T, const double mu, const double zeta_F, const double zeta_F_0){
    return A_factor(flavor, x, b_T) * B_factor(b_T, mu, zeta_F) * C_factor(flavor, x, b_T, zeta_F, zeta_F_0);
  }

  double F_kspace_integrand(const double Atan_b_T, void * par){
    double * var = (double *) par;
    int flavor = (int) (var[0] + 1.0e-9);
    double x = var[1];
    double mu = var[2];
    double zeta_F = var[3];
    double zeta_F_0 = var[4];
    double k_T = var[5];
    double b_T = tan(Atan_b_T);
    double Fb = F_bspace(flavor, x, b_T, mu, zeta_F, zeta_F_0);
    return 1.0 / (2.0 * Pi) * ROOT::Math::cyl_bessel_j(0, b_T * k_T) * b_T * Fb / pow(cos(Atan_b_T), 2);
  }

  double F_kspace(const int flavor, const double x, const double k_T, const double mu, const double zeta_F, const double zeta_F_0){
    double par[6] = {(double) flavor, x, mu, zeta_F, zeta_F_0, k_T};
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-6);
    ig.SetFunction(&F_kspace_integrand, par);
    double result = ig.Integral(0.0, Pi / 2.0 - 1.0e-9);
    return result;
  }

		    

  int Initialize(){
    order_C = 0;
    order_gamma_K = 2;
    order_gamma_F = 1;
    bstar = & bstar_CS;
    C_factor = & C_factor_AR;
    return 0;
  }
}





#endif
