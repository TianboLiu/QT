#ifndef _TMDEVOL_H_
#define _TMDEVOL_H_

#include <iostream>
#include <fstream>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

#include "TMath.h"
#include "Math/GSLIntegrator.h"

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
      if (mu > 1.3) Nf = 4.0;
      if (mu > 4.5) Nf = 5.0;
      if (mu > 180.0) Nf = 6.0;
      value += 2.0 * pow(xpdf->alphasQ(mu) / Pi, 2) * C_F / 2.0 * (C_A * (67.0 / 18.0 - Pi * Pi / 6.0) - 10.0 / 9.0 * T_R * Nf);
    }
    return value;
  }

  double gamma_F(const double mu, const double zeta_F){
    double value = (xpdf->alphasQ(mu)) * C_F / Pi * (3.0 / 2.0 - log(zeta_F / (mu * mu)));
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

  

  int Initialize(){
    bstar = & bstar_CS;
    order_C = 0;
    order_gamma_K = 1;
    order_gamma_F = 1;
    return 0;
  }
}





#endif
