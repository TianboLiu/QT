#ifndef _LCORE_H_
#define _LCORE_H_

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

#include "TString.h"


using namespace std;

const LHAPDF::PDF * pdfs = LHAPDF::mkPDF("CJ15lo", 0);
const LHAPDF::PDF * ffpion = LHAPDF::mkPDF("DSSFFlo", 211);

const double Mp = 0.938272;
const double Mpion = 0.13957;

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

  double Model_FUUT_0(const double * var, const double * par, const char * target = "proton", const char * hadron = "pi+"){
    //var: x, Q2, z, Pt
    double Pt2 = var[2] * var[2] * par[0] + par[1];
    double factor = exp(- var[3] * var[3] / Pt2) / (M_PI * Pt2);
    double xf[13], zD[13];
    DIS::xPDF(xf, var[0], var[1], target);
    DIS::zFF(zD, var[2], var[1], hadron);
    double result = factor / var[2] * (pow(2.0/3.0, 2) * (xf[2] * zD[2] + xf[4] * zD[4] + xf[6] * zD[6] + xf[2+6] * zD[2+6] + xf[4+6] * zD[4+6] + xf[6+6] * zD[6+6])
			      + pow(-1.0/3.0, 2) * (xf[1] * zD[1] + xf[3] * zD[3] + xf[5] * zD[5] + xf[1+6] * zD[1+6] + xf[3+6] * zD[3+6] + xf[5+6] * zD[5+6]));
    return result;
  }

  double Multiplicity(const double * var, const double * par, const char * target = "proton", const char * hadron = "pi+"){
    //var: x, Q2, z, Pt
    double result = 2.0 * M_PI * var[3] * Model_FUUT_0(var, par, target, hadron) / DIS::FT(var, target);
    return result;
  }

}

namespace FIT{

  int Npt = 0;
  double Value[3000], Variable[3000][4], Error[3000][2];

  double Parameters[20];

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
    TString Target[3000], Hadron[3000];
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



}











#endif
