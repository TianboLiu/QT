#ifndef _LCORE_H_
#define _LCORE_H_

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

#include "TString.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

using namespace std;

const LHAPDF::PDF * pdfs = LHAPDF::mkPDF("CJ15lo", 0);
const LHAPDF::PDF * ffpion = LHAPDF::mkPDF("DSSFFlo", 211);

const double Mp = 0.938272;
const double Mpion = 0.13957;
const double MZ = 91.1876;
const double MW = 80.385;
const double Br_Z_llbar = 0.03366;

const double alpha_EM_0 = 1.0 / 137.0;
const double alpha_EM_Z = 1.0 / 128.0;
const double sinTW2 = 0.2313;//sin theta_W square

namespace DIS{
  int Exchange(double & a, double & b){
    double tmp = b;
    b = a;
    a = tmp;
    return 0;
  }

  int xPDF(double xf[], const double x, const double Q2, const char * target = "proton"){
    xf[0] = pdfs->xfxQ2(21, x, Q2);//gluon
    for (int i = 1; i < 7; i++){
      xf[i] = pdfs->xfxQ2(i, x, Q2);//quark: d, u, s, c, b, t
      xf[i+6] = pdfs->xfxQ2(-i, x, Q2);//antiquark: db, ub, sb, cb, bb, tb
    }
    if (strcmp(target, "proton") == 0) return 0;
    else if (strcmp(target, "neutron") == 0){
      Exchange(xf[1], xf[2]);
      Exchange(xf[1+6], xf[2+6]);
      return 0;
    }
    else if (strcmp(target, "deuteron") == 0){
      double xfp[13], xfn[13];
      xPDF(xfp, x, Q2, "proton");
      xPDF(xfn, x, Q2, "neutron");
      for (int i = 0; i < 13; i++)
	xf[i] = xfp[i] + xfn[i];
      return 0;
    }
    else if (strcmp(target, "helium-3") == 0){
      double xfp[13], xfn[13];
      xPDF(xfp, x, Q2, "proton");
      xPDF(xfn, x, Q2, "neutron");
      for (int i = 0; i < 13; i++)
	xf[i] = xfp[i] * 2.0 + xfn[i];
      return 0;
    }
    else if (strcmp(target, "antiproton") == 0){
      for (int i = 1; i < 7; i++)
	Exchange(xf[i], xf[i+6]);
      return 0;
    }
    else {
      cout << "DIS::xPDF -- Wrong target option!" << endl;
      return 1;
    }
  }
      
  int zFF(double zD[], const double z, const double Q2, const char * hadron = "pi+"){
    zD[0] = ffpion->xfxQ2(21, z, Q2);
    for (int i = 0; i < 7; i++){
      zD[i] = ffpion->xfxQ2(i, z, Q2);
      zD[i+6] = ffpion->xfxQ2(-i, z, Q2);
    }
    if (strcmp(hadron, "pi+") == 0) return 0;
    else if (strcmp(hadron, "pi-") == 0){
      for (int i = 1; i < 7; i++)
	Exchange(zD[i], zD[i+6]);
      return 0;
    }
    else if (strcmp(hadron, "pi0") == 0){
      for (int i = 1; i < 7; i++){
	zD[i] = 0.5 * (zD[i] + zD[i+6]);
	zD[i+6] = zD[i];
      }
      return 0;
    }
    else {
      cout << "DIS::zFF -- Wrong hadron option!" << endl;
      return 1;
    }
  }

  double FT(const double * var, const char * target = "proton"){
    //var: x, Q2
    double xf[13];
    xPDF(xf, var[0], var[1], target);
    double result = pow(2.0/3.0, 2) * (xf[2] + xf[4] + xf[6] + xf[2+6] + xf[4+6] + xf[6+6])
      + pow(-1.0/3.0, 2) * (xf[1] + xf[3] + xf[5] + xf[1+6] + xf[3+6] + xf[5+6]);
    return result;
  }


}

namespace SIDIS{

  double (* FUUT)(const double * var, const double * par, const char * target, const char * hadron);

  double Model_FUUT_0(const double * var, const double * par, const char * target = "proton", const char * hadron = "pi+"){//simple gaussian
    //var: x, Q2, z, Pt
    double Pt2 = var[2] * var[2] * par[0] * par[0] + par[1] * par[1];
    double factor = exp(- var[3] * var[3] / Pt2) / (M_PI * Pt2);
    double xf[13], zD[13];
    DIS::xPDF(xf, var[0], var[1], target);
    DIS::zFF(zD, var[2], var[1], hadron);
    double result = factor / var[2] * (pow(2.0/3.0, 2) * (xf[2] * zD[2] + xf[4] * zD[4] + xf[6] * zD[6] + xf[2+6] * zD[2+6] + xf[4+6] * zD[4+6] + xf[6+6] * zD[6+6])
				       + pow(-1.0/3.0, 2) * (xf[1] * zD[1] + xf[3] * zD[3] + xf[5] * zD[5] + xf[1+6] * zD[1+6] + xf[3+6] * zD[3+6] + xf[5+6] * zD[5+6]));
    return result;
  }

  double Model_FUUT_1(const double * var, const double * par, const char * target = "proton", const char * hadron = "pi+"){//gaussian z-dep
    //var: x, Q2, z, Pt
    double Pt2 = var[2] * var[2] * par[0] * par[0] + par[1] * par[1] * pow(1.0 - var[2], par[2]);
    double factor = exp(- var[3] * var[3] / Pt2) / (M_PI * Pt2);                                                                                          
    double xf[13], zD[13];        
    DIS::xPDF(xf, var[0], var[1], target);
    DIS::zFF(zD, var[2], var[1], hadron);
    double result = factor / var[2] * (pow(2.0/3.0, 2) * (xf[2] * zD[2] + xf[4] * zD[4] + xf[6] * zD[6] + xf[2+6] * zD[2+6] + xf[4+6] * zD[4+6] + xf[6+6] * zD[6+6])                                                              
				       + pow(-1.0/3.0, 2) * (xf[1] * zD[1] + xf[3] * zD[3] + xf[5] * zD[5] + xf[1+6] * zD[1+6] + xf[3+6] * zD[3+6] + xf[5+6] * zD[5+6]));
    return result;
  }

  double Model_FUUT_2(const double * var, const double * par, const char * target = "proton", const char * hadron = "pi+"){//gaussian flavor-dep
    //var: x, Q2, z, Pt
    double ktv = par[0];
    double kts = par[1];
    double ptv = par[2];
    double pts = par[3];
    double z = var[2];
    double Pt2[13] = {0.0, z*z*ktv*ktv + pts*pts, z*z*ktv*ktv + ptv*ptv, z*z*kts*kts + pts*pts, z*z*kts*kts + pts*pts, z*z*kts*kts + pts*pts, z*z*kts*kts + pts*pts,
		      z*z*kts*kts + ptv*ptv, z*z*kts*kts + pts*pts, z*z*kts*kts + pts*pts, z*z*kts*kts + pts*pts, z*z*kts*kts + pts*pts, z*z*kts*kts + pts*pts};
    if (strcmp(hadron, "pi-") == 0){
      DIS::Exchange(Pt2[1], Pt2[2]);
      DIS::Exchange(Pt2[1+6], Pt2[2+6]);
    }      
    double xf[13], zD[13];
    DIS::xPDF(xf, var[0], var[1], target);
    DIS::zFF(zD, var[2], var[1], hadron);
    double factor[13];
    for (int i = 0; i < 13; i++)
      factor[i] = exp(-var[3] * var[3] / Pt2[i]) / (M_PI * Pt2[i]);
    double result = pow(2.0/3.0, 2) * (xf[2] * zD[2] * factor[2] + xf[4] * zD[4] * factor[4] + xf[6] * zD[6] * factor[6] + xf[2+6] * zD[2+6] * factor[2+6] + xf[4+6] * zD[4+6] * factor[4+6] + xf[6+6] * zD[6+6] * factor[6+6])
      + pow(-1.0/3.0, 2) * (xf[1] * zD[1] * factor[1] + xf[3] * zD[3] * factor[3] + xf[5] * zD[5] * factor[5] + xf[1+6] * zD[1+6] * factor[1+6] + xf[3+6] * zD[3+6] * factor[3+6] + xf[5+6] * zD[5+6] * factor[5+6]);
    return result / z;
  }
  
  double Multiplicity(const double * var, const double * par, const char * target = "proton", const char * hadron = "pi+"){
    //var: x, Q2, z, Pt
    //double factor = (1.0 + 2.0 * var[0] * Mp * Mp / var[1]);
    double factor = 1.0;
    double result = factor * 2.0 * M_PI * var[3] * FUUT(var, par, target, hadron) / DIS::FT(var, target);
    return result;
  }

}

namespace DY{
  
  double (* FUU1DY)(const double * var, const double * par, const char * hadron1, const char * hadron2);

  double Model_FUU1DY_0(const double * var, const double * par, const char * hadron1 = "proton", const char * hadron2 = "proton"){
    //var: Q2, qT, y, s
    double x1 = sqrt(var[0] / var[3]) * exp(var[2]);
    double x2 = sqrt(var[0] / var[3]) * exp(-var[2]);
    double xf1[13], xf2[13];
    DIS::xPDF(xf1, x1, var[0], hadron1);
    DIS::xPDF(xf2, x2, var[0], hadron2);
    double factor = (pow(2.0 / 3.0, 2) * (xf1[2] * xf2[2+6] + xf1[4] * xf2[4+6] + xf1[6] * xf2[6+6] + xf1[2+6] * xf2[2] + xf1[4+6] * xf2[4] + xf1[6+6] * xf2[6]) + pow(1.0 / 3.0, 2) * (xf1[1] * xf2[1+6] + xf1[3] * xf2[3+6] + xf1[5] * xf2[5+6] + xf1[1+6] * xf2[1] + xf1[3+6] * xf2[3] + xf1[5+6] * xf2[5])) / (x1 * x2);
    double kt2 = par[0] * par[0];
    return factor / 3.0 * exp(-var[1] * var[1] / (2.0 * kt2)) / (2.0 * M_PI * kt2);
  }

  double (* FUU1Z)(const double * var, const double * par, const char * hadron1, const char * hadron2);

  double Model_FUU1Z_0(const double * var, const double * par, const char * hadron1 = "proton", const char * hadron2 = "proton"){
    //var: Q2, qT, y, s
    double x1 = sqrt(var[0] / var[3]) * exp(var[2]);
    double x2 = sqrt(var[0] / var[3]) * exp(-var[2]);
    double xf1[13], xf2[13];
    DIS::xPDF(xf1, x1, var[0], hadron1);
    DIS::xPDF(xf2, x2, var[0], hadron2);
    double factor = ((pow(0.5 - 2.0 * 2.0 / 3.0 * sinTW2, 2) + pow(0.5, 2)) * (xf1[2] * xf2[2+6] + xf1[4] * xf2[4+6] + xf1[6] * xf2[6+6] + xf1[2+6] * xf2[2] + xf1[4+6] * xf2[4] + xf1[6+6] * xf2[6])
		     + (pow(0.5 + 2.0 * 1.0 / 3.0 * sinTW2, 2) + pow(-0.5, 2)) * (xf1[1] * xf2[1+6] + xf1[3] * xf2[3+6] + xf1[5] * xf2[5+6] + xf1[1+6] * xf2[1] + xf1[3+6] * xf2[3] + xf1[5+6] * xf2[5])) / (x1 * x2);
    double kt2 = par[0] * par[0];
    return factor / 3.0 * exp(-var[1] * var[1] / (2.0 * kt2)) / (2.0 * M_PI * kt2);
  }
  
}



namespace FIT{

  int Npt = 0;
  double Value[3000], Variable[3000][4], Error[3000][2];
  TString Target[3000], Hadron[3000];

  double Parameters[20];
  double ParametersError[20];

  double SelectionT[4];
  double SelectionTdelta[4];

  int PrintLevel = 1;

  int CheckValue(const double x[], const double t[], const double dt[]){
    if (pow(x[0] - t[0], 2) < pow(dt[0], 2) &&
	pow(x[1] - t[1], 2) < pow(dt[1], 2) && 
	pow(x[2] - t[2], 2) < pow(dt[2], 2) && 
	pow(x[3] - t[3], 2) < pow(dt[3], 2)
	) return 1;
    else return 0;
  }

  int LoadData(const char * filename, const int skiprows = 19, const char * target = "proton", const char * hadron = "pi+"){
    ifstream infile(filename);
    char ltmp[300];
    for (int i = 0; i < skiprows; i++)
      infile.getline(ltmp, 300);
    double var[4], value, error[2], tmp;
    while (infile >> tmp >> value >> error[0] >> error[1] >> var[1] >> var[0] >> var[2] >> var[3]){
      if (CheckValue(var, SelectionT, SelectionTdelta)){
	for (int i = 0; i < 4; i++)
	  Variable[Npt][i] = var[i];
	Value[Npt] = value;
	Error[Npt][0] = error[0];
	Error[Npt][1] = error[1];
	Target[Npt] = target;
	Hadron[Npt] = hadron;
	if (PrintLevel > 0)
	  cout << "  " << Npt << " " << var[0] << " " << var[1] << " " << var[2] << " " << var[3] << " " << value << " " << error[0] << " " << error[1] << " " << Target[Npt].Data() << " " << Hadron[Npt].Data() << endl;
	Npt++;
      }
    }
    infile.close();
    return 0;
  }
  
  double Chi2(const double * par){
    double sum = 0.0;
    for (int i = 0; i < Npt; i++){
      sum += pow(SIDIS::Multiplicity(Variable[i], par, Target[i].Data(), Hadron[i].Data()) - Value[i], 2) / (pow(Error[i][0], 2) + pow(Error[i][1], 2));
    }
    return sum;
  }

  double Minimize(const int NPAR, const double * init){
    ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    min->SetMaxFunctionCalls(100000);
    //min->SetPrecision(1.0e-14);
    min->SetTolerance(1.0e-3);
    min->SetPrintLevel(0);
    ROOT::Math::Functor f(&Chi2, NPAR);
    min->SetFunction(f);
    for (int i = 0; i < NPAR; i++){
      min->SetVariable(i, Form("p%d", i), init[i], 1.0e-4);
    }
    min->Minimize();
    const double * xs = min->X();
    const double * es = min->Errors();
    for (int i = 0; i < NPAR; i++){
      Parameters[i] = xs[i];
      ParametersError[i] = es[i];
    }
    const double chi2 = min->MinValue();
    Parameters[NPAR] = chi2;
    return chi2;
  }

  int Compare(){
    for (int i = 0; i < Npt; i++){
      cout << SIDIS::Multiplicity(Variable[i], Parameters, Target[i].Data(), Hadron[i].Data()) << " -- " << Value[i] << endl;
    }
    return 0;
  }


}

namespace FITDY{

  int Npt = 0;
  double Value[3000], Variable[3000][4], Error[3000][2];
  int Label[3000];

  double Parameters[20];
  double ParametersError[20];

  double SelectionT[4];
  double SelectionTdelta[4];

  int PrintLevel = 1;

  int CheckValue(const double x[], const double t[], const double dt[]){
    if (pow(x[0] - t[0], 2) < pow(dt[0], 2) &&
	pow(x[1] - t[1], 2) < pow(dt[1], 2) && 
	pow(x[2] - t[2], 2) < pow(dt[2], 2) && 
	pow(x[3] - t[3], 2) < pow(dt[3], 2)
	) return 1;
    else return 0;
  }
  
  int LoadData(const char * filename, const int skiprows = 19, const char * experiment = "E288_200"){
    ifstream infile(filename);
    char ltmp[300];
    for (int i = 0; i < skiprows; i++)
      infile.getline(ltmp, 300);
    int label = 0;
    double s = 0.0;
    if (strcmp(experiment, "E288_200") == 0){
      //cout << "Loading E288 200 data" << endl;
      label = 1;
      s = pow(200.0 + Mp, 2) - pow(200.0, 2);
      double var[4], value, error[2];
      var[2] = 0.40;
      var[3] = s;
      while (infile >> var[0] >> var[1] >> value >> error[0]){
	if (CheckValue(var, SelectionT, SelectionTdelta)){
	  Variable[Npt][0] = var[0] * var[0];
	  Variable[Npt][1] = var[1];
	  Variable[Npt][2] = var[2];
	  Variable[Npt][3] = var[3];
	  Value[Npt] = value * pow(1.0e13 / 0.197, 2);
	  Error[Npt][0] = error[0] * pow(1.0e13 / 0.197, 2);
	  Error[Npt][1] = 0.25 * Value[Npt];
	  Label[Npt] = label;
	  if (PrintLevel > 0){
	    cout << "  " << Npt << " " << Variable[Npt][0] << " " << Variable[Npt][1] << " " << Variable[Npt][2] << " " << Variable[Npt][3] << " " << Value[Npt] << " " << Error[Npt][0] << Error[Npt][1] << endl;
	  }
	  Npt++;
	}
      }
      infile.close();
      return 0;
    }
    if (strcmp(experiment, "E288_300") == 0){
      //cout << "Loading E288 300 data" << endl;
      label = 1;
      s = pow(300.0 + Mp, 2) - pow(300.0, 2);
      double var[4], value, error[2];
      var[2] = 0.21;
      var[3] = s;
      while (infile >> var[0] >> var[1] >> value >> error[0]){
	if (CheckValue(var, SelectionT, SelectionTdelta)){
	  Variable[Npt][0] = var[0] * var[0];
	  Variable[Npt][1] = var[1];
	  Variable[Npt][2] = var[2];
	  Variable[Npt][3] = var[3];
	  Value[Npt] = value * pow(1.0e13 / 0.197, 2);
	  Error[Npt][0] = error[0] * pow(1.0e13 / 0.197, 2);
	  Error[Npt][1] = 0.25 * Value[Npt];
	  Label[Npt] = label;
	  if (PrintLevel > 0){
	    cout << "  " << Npt << " " << Variable[Npt][0] << " " << Variable[Npt][1] << " " << Variable[Npt][2] << " " << Variable[Npt][3] << " " << Value[Npt] << " " << Error[Npt][0] << Error[Npt][1] << endl;
	  }
	  Npt++;
	}
      }
      infile.close();
      return 0;
    }
    if (strcmp(experiment, "E288_400") == 0){
      //cout << "Loading E288 400 data" << endl;
      label = 1;
      s = pow(400.0 + Mp, 2) - pow(400.0, 2);
      double var[4], value, error[2];
      var[2] = 0.03;
      var[3] = s;
      while (infile >> var[0] >> var[1] >> value >> error[0]){
	if (CheckValue(var, SelectionT, SelectionTdelta)){
	  Variable[Npt][0] = var[0] * var[0];
	  Variable[Npt][1] = var[1];
	  Variable[Npt][2] = var[2];
	  Variable[Npt][3] = var[3];
	  Value[Npt] = value * pow(1.0e13 / 0.197, 2);
	  Error[Npt][0] = error[0] * pow(1.0e13 / 0.197, 2);
	  Error[Npt][1] = 0.25 * Value[Npt];
	  Label[Npt] = label;
	  if (PrintLevel > 0){
	    cout << "  " << Npt << " " << Variable[Npt][0] << " " << Variable[Npt][1] << " " << Variable[Npt][2] << " " << Variable[Npt][3] << " " << Value[Npt] << " " << Error[Npt][0] << Error[Npt][1] << endl;
	  }
	  Npt++;
	}
      }
      infile.close();
      return 0;
    }
    if (strcmp(experiment, "E605_800") == 0){
      //cout << "Loading E605 800 data" << endl;
      label = 2;
      s = pow(800.0 + Mp, 2) - pow(800.0, 2);
      double var[4], value, error[2];
      var[2] = 0.0;
      var[3] = s;
      while (infile >> var[0] >> var[1] >> value >> error[0]){
	var[2] = asinh(sqrt(s) / var[0] * 0.1 / 2.0);//x_F = 0.1
	if (CheckValue(var, SelectionT, SelectionTdelta)){
	  Variable[Npt][0] = var[0] * var[0];
	  Variable[Npt][1] = var[1];
	  Variable[Npt][2] = var[2];
	  Variable[Npt][3] = var[3];
	  Value[Npt] = value * pow(1.0e13 / 0.197, 2);
	  Error[Npt][0] = error[0] * pow(1.0e13 / 0.197, 2);
	  Error[Npt][1] = 0.15 * Value[Npt];
	  Label[Npt] = label;
	  if (PrintLevel > 0){
	    cout << "  " << Npt << " " << Variable[Npt][0] << " " << Variable[Npt][1] << " " << Variable[Npt][2] << " " << Variable[Npt][3] << " " << Value[Npt] << " " << Error[Npt][0] << Error[Npt][1] << endl;
	  }
	  Npt++;
	}
      }
      infile.close();
      return 0;
    }
    return 0;
  }
  
  double Chi2(const double * par){
    double sum = 0.0;
    double theory = 0.0;
    for (int i = 0; i < Npt; i++){
      if (Label[i] == 1){
	theory = (4.0 * M_PI * pow(alpha_EM_0, 2) / (3.0 * Variable[i][3]) * 2.0 * log((sqrt(Variable[i][0]) + 0.5) / (sqrt(Variable[i][0]) - 0.5)))
	  * (78.0/195.0 * DY::FUU1DY(Variable[i], par, "proton", "proton") + 117.0/195.0 * DY::FUU1DY(Variable[i], par, "proton", "neutron"));
      }
      else if (Label[i] == 2){
	if (Variable[i][0] < pow(11.5, 2)){
	  theory = (4.0 * M_PI * pow(alpha_EM_0, 2) / (3.0 * Variable[i][3]) * 2.0 * log((sqrt(Variable[i][0]) + 0.5) / (sqrt(Variable[i][0]) - 0.5)))
	    * (29.0/63.0 * DY::FUU1DY(Variable[i], par, "proton", "proton") + 34.0/63.0 * DY::FUU1DY(Variable[i], par, "proton", "neutron"));
	}
	else if (Variable[i][0] < pow(13.5, 2)){
	  theory = (4.0 * M_PI * pow(alpha_EM_0, 2) / (3.0 * Variable[i][3]) * 2.0 * log((sqrt(Variable[i][0]) + 1.0) / (sqrt(Variable[i][0]) - 1.0)))
	    * (29.0/63.0 * DY::FUU1DY(Variable[i], par, "proton", "proton") + 34.0/63.0 * DY::FUU1DY(Variable[i], par, "proton", "neutron"));
	}
	else {
	  theory = (4.0 * M_PI * pow(alpha_EM_0, 2) / (3.0 * Variable[i][3]) * 2.0 * log((sqrt(Variable[i][0]) + 2.2) / (sqrt(Variable[i][0]) - 2.3)))
	    * (29.0/63.0 * DY::FUU1DY(Variable[i], par, "proton", "proton") + 34.0/63.0 * DY::FUU1DY(Variable[i], par, "proton", "neutron"));
	}
      } 
      else {
	cout << "Wrong Label" << endl;
	return -1;
      }
      sum += pow(theory - Value[i], 2) / (pow(Error[i][0], 2) + pow(Error[i][1], 2));
    }
    return sum;
  }

  double Minimize(const int NPAR, const double * init){
    ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    min->SetMaxFunctionCalls(100000);
    //min->SetPrecision(1.0e-14);
    min->SetTolerance(1.0e-3);
    min->SetPrintLevel(0);
    ROOT::Math::Functor f(&Chi2, NPAR);
    min->SetFunction(f);
    for (int i = 0; i < NPAR; i++){
      min->SetVariable(i, Form("p%d", i), init[i], 1.0e-4);
    }
    min->Minimize();
    const double * xs = min->X();
    const double * es = min->Errors();
    for (int i = 0; i < NPAR; i++){
      Parameters[i] = xs[i];
      ParametersError[i] = es[i];
    }
    const double chi2 = min->MinValue();
    Parameters[NPAR] = chi2;
    return chi2;
  }

 
}

namespace FITZ{

  int Npt = 0;
  double Value[3000], Variable[3000][4], Error[3000][2];
  int Label[3000];

  double Parameters[20];
  double ParametersError[20];

  double SelectionT[4];
  double SelectionTdelta[4];

  int PrintLevel = 1;

  int CheckValue(const double x[], const double t[], const double dt[]){
    if (pow(x[0] - t[0], 2) < pow(dt[0], 2) &&
	pow(x[1] - t[1], 2) < pow(dt[1], 2) && 
	pow(x[2] - t[2], 2) < pow(dt[2], 2) && 
	pow(x[3] - t[3], 2) < pow(dt[3], 2)
	) return 1;
    else return 0;
  }
  
  int LoadData(const char * filename, const int skiprows = 19, const char * experiment = "CDF_RunI"){
    ifstream infile(filename);
    char ltmp[300];
    for (int i = 0; i < skiprows; i++)
      infile.getline(ltmp, 300);
    int label = 0;
    double s = 0.0;
    if (strcmp(experiment, "CDF_RunI") == 0){
      //cout << "Loading CDF RunI data" << endl;
      label = 3;
      s = pow(1.8e3, 2);
      double var[4], value, error[2];
      var[0] = MZ;
      var[3] = s;
      while (infile >> var[1] >> value >> error[0]){
	if (CheckValue(var, SelectionT, SelectionTdelta)){
	  Variable[Npt][0] = var[0] * var[0];
	  Variable[Npt][1] = var[1];
	  Variable[Npt][3] = var[3];
	  Value[Npt] = value * 1.0e-10 * pow(1.0 / 0.197, 2);
	  Error[Npt][0] = error[0] * 1.0e-10 * pow(1.0 / 0.197, 2);
	  Error[Npt][1] = 0.0;
	  Label[Npt] = label;
	  if (PrintLevel > 0){
	    cout << "  " << Npt << " " << Variable[Npt][0] << " " << Variable[Npt][1] << " " << Variable[Npt][3] << " " << Value[Npt] << " " << Error[Npt][0] << Error[Npt][1] << endl;
	  }
	  Npt++;
	}
      }
      infile.close();
      return 0;
    }
    if (strcmp(experiment, "CDF_RunII") == 0){
      //cout << "Loading CDF RunII data" << endl;
      label = 3;
      s = pow(1.96e3, 2);
      double var[4], value, error[2];
      var[0] = MZ;
      var[3] = s;
      while (infile >> var[1] >> value >> error[0] >> error[1]){
	if (CheckValue(var, SelectionT, SelectionTdelta)){
	  Variable[Npt][0] = var[0] * var[0];
	  Variable[Npt][1] = var[1];
	  Variable[Npt][3] = var[3];
	  Value[Npt] = value * 1.0e-10 * pow(1.0 / 0.197, 2);
	  Error[Npt][0] = error[0] * 1.0e-10 * pow(1.0 / 0.197, 2);
	  Error[Npt][1] = error[1] * 1.0e-10 * pow(1.0 / 0.197, 2);
	  Label[Npt] = label;
	  if (PrintLevel > 0){
	    cout << "  " << Npt << " " << Variable[Npt][0] << " " << Variable[Npt][1] << " " << Variable[Npt][3] << " " << Value[Npt] << " " << Error[Npt][0] << Error[Npt][1] << endl;
	  }
	  Npt++;
	}
      }
      infile.close();
      return 0;
    }
    if (strcmp(experiment, "D0_RunI") == 0){
      //cout << "Loading D0 RunI data" << endl;
      label = 3;
      s = pow(1.8e3, 2);
      double var[4], value, error[2];
      var[0] = MZ;
      var[3] = s;
      while (infile >> var[1] >> value >> error[0]){
	if (CheckValue(var, SelectionT, SelectionTdelta)){
	  Variable[Npt][0] = var[0] * var[0];
	  Variable[Npt][1] = var[1];
	  Variable[Npt][3] = var[3];
	  Value[Npt] = value * 1.0e-10 * pow(1.0 / 0.197, 2);
	  Error[Npt][0] = error[0] * 1.0e-10 * pow(1.0 / 0.197, 2);
	  Error[Npt][1] = 0.0;
	  Label[Npt] = label;
	  if (PrintLevel > 0){
	    cout << "  " << Npt << " " << Variable[Npt][0] << " " << Variable[Npt][1] << " " << Variable[Npt][3] << " " << Value[Npt] << " " << Error[Npt][0] << Error[Npt][1] << endl;
	  }
	  Npt++;
	}
      }
      infile.close();
      return 0;
    }
    if (strcmp(experiment, "D0_RunII") == 0){
      //cout << "Loading D0 RunII data" << endl;
      label = 3;
      s = pow(1.96e3, 2);
      double var[4], value, error[2];
      var[0] = MZ;
      var[3] = s;
      while (infile >> var[1] >> value >> error[0] >> error[1]){
	if (CheckValue(var, SelectionT, SelectionTdelta)){
	  Variable[Npt][0] = var[0] * var[0];
	  Variable[Npt][1] = var[1];
	  Variable[Npt][3] = var[3];
	  Value[Npt] = 255.8 * value * 1.0e-10 * pow(1.0 / 0.197, 2);
	  Error[Npt][0] = Value[Npt] * sqrt((pow(error[0], 2) + pow(error[1], 2)) / pow(value, 2));
	  Error[Npt][1] = Value[Npt] * 16.0 / 255.8;
	  Label[Npt] = label;
	  if (PrintLevel > 0){
	    cout << "  " << Npt << " " << Variable[Npt][0] << " " << Variable[Npt][1] << " " << Variable[Npt][3] << " " << Value[Npt] << " " << Error[Npt][0] << Error[Npt][1] << endl;
	  }
	  Npt++;
	}
      }
      infile.close();
      return 0;
    }
    return 0;
  }

  double IntegrandPPbar(const double y, void * others){
    double * var = (double *) others;
    double variable[4] = {var[0], var[1], y, var[3]};
    double result = DY::FUU1Z(variable, &var[4], "proton", "antiproton");
    return result;
  }
  
  double Chi2(const double * par){
    double sum = 0.0;
    double theory = 0.0;
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-6);
    double others[24];
    ig.SetFunction(&S_integrand, &others);
    for (int i = 0; i < Npt; i++){
      if (Label[i] == 3){
	others[0] = Variable[i][0];
	others[1] = Variable[i][1];
	others[3] = Variable[i][3];
	for (int j = 0; j < 20; j++)
	  others[j+4] = par[j];
	theory = ig.Integral(0.5 * log(Variable[i][0] / Variable[i][3]), 0.5 * log(Variable[i][3] / Variable[i][0]))
	  * 2.0 * Variable[i][1] * (M_PI * M_PI * alpha_EM_Z / (Variable[i][3] * sinTW2 * (1.0 - sinTW2))) * 0.03366;
      }
      else {
	cout << "Wrong Label" << endl;
	return -1;
      }
      sum += pow(theory - Value[i], 2) / (pow(Error[i][0], 2) + pow(Error[i][1], 2));
    }
    return sum;
  }

  double Minimize(const int NPAR, const double * init){
    ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    min->SetMaxFunctionCalls(100000);
    //min->SetPrecision(1.0e-14);
    min->SetTolerance(1.0e-3);
    min->SetPrintLevel(0);
    ROOT::Math::Functor f(&Chi2, NPAR);
    min->SetFunction(f);
    for (int i = 0; i < NPAR; i++){
      min->SetVariable(i, Form("p%d", i), init[i], 1.0e-4);
    }
    min->Minimize();
    const double * xs = min->X();
    const double * es = min->Errors();
    for (int i = 0; i < NPAR; i++){
      Parameters[i] = xs[i];
      ParametersError[i] = es[i];
    }
    const double chi2 = min->MinValue();
    Parameters[NPAR] = chi2;
    return chi2;
  }

 
}








#endif
