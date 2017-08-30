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
    double value = 2.0 * (xpdf->alphasQ(mu)) * C_F / Pi;
    return value;
  }

  double gamma_F(const double mu, const double zeta_F){
    double value = (xpdf->alphasQ(mu)) * C_F / Pi * (3.0 / 2.0 - log(zeta_F / (mu * mu)));
    return value;
  }

  double logB_integrand_2(const double logmu, const double * par){
    double zeta_F = par[0];
    double mu = exp(logmu);
    return gamma_F(mu, mu * mu) - log(sqrt(zeta_F) / mu) * gamma_K(mu);
  }

  double logB_2(const double b_T, const double mu, const double zeta_F){
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-6);
    ig.SetFunction(&logB_integrand_2, &zeta_F);
    double min_logmu = log(mu_b(b_T));
    double max_logmu = log(mu);
    double result = ig.Integral(min_logmu, max_logmu);
    return result;
  }


  int Initialize(){
    bstar = & bstar_CS;
    return 0;
  }
}





#endif
