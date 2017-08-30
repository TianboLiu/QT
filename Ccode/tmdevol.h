#ifndef _TMDEVOL_H_
#define _TMDEVOL_H_

#include <iostream>
#include <fstream>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

#include "TMath.h"

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
  

}

TMDEVOL::bstar = & TMDEVOL::bstar_CS;



#endif
