#ifndef _THEORY_H_
#define _THEORY_H_

#include <iostream>
#include <fstream>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

#include "Math/GSLIntegrator.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Math/WrappedTF1.h"
#include "Math/WrappedParamFunction.h"
#include "TF1.h"
#include "Math/SpecFuncMathMore.h"
#include "gsl/gsl_sf_bessel.h"

using namespace std;

double J0(const double x){//Regular Bessel funtion
  return ROOT::Math::cyl_bessel_j(0, x)
  //return gsl_sf_bessel_J0(x);
} 

const LHAPDF::PDF * xpdf = LHAPDF::mkPDF("CJ15lo", 0);

double (* F_TMD)(const int flavor, const double x, const double b_T, const double Q);

double Parameters[3] = {1.0, 1.0, 1.0};
double F_TMD_simple(const int flavor, const double x, const double b_T, const double Q){
  double col = xpdf->xfxQ(flavor, x, Q) / x;
  double trans = exp(-Parameters[2] * pow(b_T, Parameters[1]));
  return Parameters[0] * col * trans;
}

int GetTMDs(double * tmds, const double dProton, const double x, const double b_T, const double Q){
  double dP = abs(dProton);
  double dN = 1.0 - dP;
  tmds[2] = dP * F_TMD(2, x, b_T, Q) + dN * F_TMD(1, x, b_T, Q);
  tmds[1] = dP * F_TMD(1, x, b_T, Q) + dN * F_TMD(2, x, b_T, Q);
  tmds[2+6] = dP * F_TMD(-2, x, b_T, Q) + dN * F_TMD(-1, x, b_T, Q);
  tmds[1+6] = dP * F_TMD(-1, x, b_T, Q) + dN * F_TMD(-2, x, b_T, Q);
  tmds[0] = F_TMD(21, x, b_T, Q);
  for (int i = 3; i <= 6; i++){
    tmds[i] = F_TMD(i, x, b_T, Q);
    tmds[i+6] = F_TMD(-i, x, b_T, Q);
  }
  if (dProton < 0){
    double tmp = 0;
    for (int i = 1; i <= 6; i++){
      tmp = tmds[i];
      tmds[i] = tmds[i+6];
      tmds[i+6] = tmp;
    }
  }
  return 0;
}

double dsigma_DY_integrand(const double * atb, const double * para){
  double b_T = tan(atb[0]);
  double Q = para[0];
  double QT = para[1];
  double y = para[2];
  double s = para[3];
  double dProtonA = para[4];
  double dProtonB = para[5];
  double xA = Q / sqrt(s) * exp(y);
  double xB = Q / sqrt(s) * exp(-y);
  if (xA > 1.0 || xB > 1.0) return 0;
  double alpha_EM_0 = 1.0 / 137.0;
  double e_u = 2.0 / 3.0;
  double e_d = -1.0 / 3.0;
  double sigma0 = 4.0 * pow(M_PI, 2) * pow(alpha_EM_0, 2) / (9.0 * s * pow(Q, 2));
  double tmdA[13], tmdB[13];
  GetTMDs(tmdA, dProtonA, xA, b_T, Q);
  GetTMDs(tmdB, dProtonB, xB, b_T, Q);
  double convol =
    pow(e_u, 2) * (tmdA[2] * tmdB[2+6] + tmdA[2+6] * tmdB[2]
		   + tmdA[4] * tmdB[4+6] + tmdA[4+6] * tmdB[4])
    +
    pow(e_d, 2) * (tmdA[1] * tmdB[1+6] + tmdA[1+6] * tmdB[1]
		   + tmdA[3] * tmdB[3+6] + tmdA[3+6] * tmdB[3]
		   + tmdA[5] * tmdB[5+6] + tmdA[5+6] * tmdB[5]);
  double result = b_T * J0(b_T * QT) * sigma0 * convol / (2.0 * M_PI);
  return result / pow(cos(atb[0]), 2);
}

double dsigma_DY(double * par){
  //par: Q, QT, y, s, dPA, dPB
  TF1 f("integrand", &dsigma_DY_integrand, 0.0, M_PI_2, 6);
  ROOT::Math::WrappedTF1 wf(f);
  wf.SetParameters(par);
  ROOT::Math::IntegratorOneDim ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1e-4);
  ig.SetFunction(wf);
  double result = ig.Integral(0.0, M_PI_2 - 0.1) + ig.Integral(M_PI_2 - 0.1, M_PI_2 - 0.01);
  return result;
}

double dsigma_Z_integrand(const double * atb, const double * para){
  double b_T = tan(atb[0]);
  double Q = para[0];
  double QT = para[1];
  double y = para[2];
  double s = para[3];
  double dProtonA = para[4];
  double dProtonB = para[5];
  double xA = Q / sqrt(s) * exp(y);
  double xB = Q / sqrt(s) * exp(-y);
  if (xA > 1.0 || xB > 1.0) return 0;
  double alpha_EM_Z = 1.0 / 128.0;
  double MZ = 91.2;
  double GammaZ = 2.5;
  double sw2 = 0.2313;
  double cw2 = 0.7687;
  double cu2 = pow(Q, 4) / (pow(Q * Q - MZ * MZ, 2) + pow(GammaZ * MZ, 2)) * (1.0 - 4.0 * sw2 + 8.0 * pow(sw2, 2)) / (8.0 * sw2 * cw2) * (1.0 - 4.0 * 2.0 / 3.0 * sw2 + 8.0 * pow(2.0 / 3.0, 2) * pow(sw2, 2)) / (8.0 * sw2 * cw2);
  double cd2 = pow(Q, 4) / (pow(Q * Q - MZ * MZ, 2) + pow(GammaZ * MZ, 2)) * (1.0 - 4.0 * sw2 + 8.0 * pow(sw2, 2)) / (8.0 * sw2 * cw2) * (1.0 - 4.0 * 1.0 / 3.0 * sw2 + 8.0 * pow(1.0 / 3.0, 2) * pow(sw2, 2)) / (8.0 * sw2 * cw2);
  double sigma0 = 4.0 * pow(M_PI, 2) * pow(alpha_EM_Z, 2) / (9.0 * s * pow(Q, 2));
  double tmdA[13], tmdB[13];
  GetTMDs(tmdA, dProtonA, xA, b_T, Q);
  GetTMDs(tmdB, dProtonB, xB, b_T, Q);
  double convol =
    cu2 * (tmdA[2] * tmdB[2+6] + tmdA[2+6] * tmdB[2]
		   + tmdA[4] * tmdB[4+6] + tmdA[4+6] * tmdB[4])
    +
    cd2 * (tmdA[1] * tmdB[1+6] + tmdA[1+6] * tmdB[1]
		   + tmdA[3] * tmdB[3+6] + tmdA[3+6] * tmdB[3]
		   + tmdA[5] * tmdB[5+6] + tmdA[5+6] * tmdB[5]);
  double result = b_T * J0(b_T * QT) * sigma0 * convol / (2.0 * M_PI);
  if (b_T * QT > 200.0) return 0;
  return result / pow(cos(atb[0]), 2);
}

double dsigma_Z(double * par){
  //par: Q, QT, y, s, dPA, dPB
  TF1 f("integrand", &dsigma_Z_integrand, 0.0, M_PI_2, 6);
  ROOT::Math::WrappedTF1 wf(f);
  wf.SetParameters(par);
  ROOT::Math::IntegratorOneDim ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1e-6);
  ig.SetFunction(wf);
  double result = ig.Integral(0.0, M_PI_2);
  return result;
}

double dsigma_gZ_integrand(const double * atb, const double * para){
  double b_T = tan(atb[0]);
  double Q = para[0];
  double QT = para[1];
  double y = para[2];
  double s = para[3];
  double dProtonA = para[4];
  double dProtonB = para[5];
  double xA = Q / sqrt(s) * exp(y);
  double xB = Q / sqrt(s) * exp(-y);
  if (xA > 1.0 || xB > 1.0) return 0;
  double alpha_EM_Z = 1.0 / 128.0;
  double MZ = 91.2;
  double GammaZ = 2.5;
  double sw2 = 0.2313;
  double cw2 = 0.7687;
  double cu2 = 2.0 * pow(Q, 2) * (Q * Q - MZ * MZ) / (pow(Q * Q - MZ * MZ, 2) + pow(GammaZ * MZ, 2)) * (1.0 - 4.0 * sw2) / (4.0 * sqrt(sw2 * cw2)) * 2.0 / 3.0 * (1.0 - 4.0 * 2.0 / 3.0 * sw2) / (4.0 * sqrt(sw2 * cw2));
  double cd2 = 2.0 * pow(Q, 2) * (Q * Q - MZ * MZ) / (pow(Q * Q - MZ * MZ, 2) + pow(GammaZ * MZ, 2)) * (1.0 - 4.0 * sw2) / (4.0 * sqrt(sw2 * cw2)) * 1.0 / 3.0 * (1.0 - 4.0 * 1.0 / 3.0 * sw2) / (4.0 * sqrt(sw2 * cw2));
  double sigma0 = 4.0 * pow(M_PI, 2) * pow(alpha_EM_Z, 2) / (9.0 * s * pow(Q, 2));
  double tmdA[13], tmdB[13];
  GetTMDs(tmdA, dProtonA, xA, b_T, Q);
  GetTMDs(tmdB, dProtonB, xB, b_T, Q);
  double convol =
    cu2 * (tmdA[2] * tmdB[2+6] + tmdA[2+6] * tmdB[2]
		   + tmdA[4] * tmdB[4+6] + tmdA[4+6] * tmdB[4])
    +
    cd2 * (tmdA[1] * tmdB[1+6] + tmdA[1+6] * tmdB[1]
		   + tmdA[3] * tmdB[3+6] + tmdA[3+6] * tmdB[3]
		   + tmdA[5] * tmdB[5+6] + tmdA[5+6] * tmdB[5]);
  double result = b_T * J0(b_T * QT) * sigma0 * convol / (2.0 * M_PI);
  if (b_T * QT > 200.0) return 0;
  return result / pow(cos(atb[0]), 2);
}

double dsigma_gZ(double * par){
  //par: Q, QT, y, s, dPA, dPB
  TF1 f("integrand", &dsigma_gZ_integrand, 0.0, M_PI_2, 6);
  ROOT::Math::WrappedTF1 wf(f);
  wf.SetParameters(par);
  ROOT::Math::IntegratorOneDim ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1e-6);
  ig.SetFunction(wf);
  double result = ig.Integral(0.0, M_PI_2);
  return result;
}


#endif
