#ifndef _LOAD_H_
#define _LOAD_H_

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

const double Mp = 0.938272;

int Npt = 0;

double Variable[2000][4], Value[2000], Errors[2000][2];
int Observable[2000];
bool FlagQ[2000];

/* Observable list
   0: DY p Pt cross section Edsigma/d3p, [Q, QT, y, s]
   1: DY p Cu cross section Edsigma/d3p, [Q, QT, y, s]
   2: DY p d  cross section Edsigma/d3p, [Q, QT, y, s]
   3:

 */

int SetFlagQ(const double Q, const double dQ){
  for (int i = 0; i < Npt; i++){
    if (pow(Variable[i][0] - Q, 2) < pow(dQ, 2)){
      FlagQ[i] = true;
    }
    else {
      FlagQ[i] = false;
    }
  }
  return 0;
}

int LoadData_DY(const char * filename, const char * experiment){
  ifstream infile(filename);
  char ltmp[300];
  int obs;
  double var[4], value, error[2];
  //E288
  if (strcmp(experiment, "E288_200") == 0 || (strcmp(experiment, "E288_300") == 0 || strcmp(experiment, "E288_400") == 0)){
    for (int i = 0; i < 15; i++)//skiprows
      infile.getline(ltmp, 300);
    obs = 0;
    if (strcmp(experiment, "E288_200") == 0){
      var[2] = 0.40;//y, rapidity
      var[3] = pow(200.0 + Mp, 2) - pow(200.0, 2);//s
    }
    else if (strcmp(experiment, "E288_300") == 0){
      var[2] = 0.21;//y, rapidity
      var[3] = pow(300.0 + Mp, 2) - pow(300.0, 2);//s
    }
    else if (strcmp(experiment, "E288_400") == 0){
      var[2] = 0.03;//y, rapidity
      var[3] = pow(400.0 + Mp, 2) - pow(400.0, 2);//s
    } 
    while (infile >> var[0] >> var[1] >> value >>  error[0]){
      Variable[Npt][0] = var[0];
      Variable[Npt][1] = var[1];
      Variable[Npt][2] = var[2];
      Variable[Npt][3] = var[3];
      Value[Npt] = value * pow(1.0e13 / 0.197327, 2);
      Errors[Npt][0] = error[0] * pow(1.0e13 / 0.197327, 2);
      Errors[Npt][1] = 0.25 * Value[Npt];
      Observable[Npt] = obs;
      Npt++;
    }
  }
  //E605
  if (strcmp(experiment, "E605") == 0){
    for (int i = 0; i < 15; i++)//skiprows
      infile.getline(ltmp, 300);
    obs = 1;
    var[3] = pow(800.0 + Mp, 2) - pow(800.0, 2);
    double xF = 0.1;
    while (infile >> var[0] >> var[1] >> value >> error[0]){
      Variable[Npt][0] = var[0];
      Variable[Npt][1] = var[1];
      Variable[Npt][2] = asinh(sqrt(var[3]) / var[0] * xF / 2.0);
      Variable[Npt][3] = var[3];
      Value[Npt] = value * pow(1.0e13 / 0.197327, 2);
      Errors[Npt][0] = error[0] * pow(1.0e13 / 0.197327, 2);
      Errors[Npt][1] = 0.15 * Value[Npt];
      Observable[Npt] = obs;
      Npt++;
    }
  }
  //E772
  if (strcmp(experiment, "E772") == 0){
    for (int i = 0; i < 15; i++)//skiprows
      infile.getline(ltmp, 300);
    obs = 2;
    var[3] = pow(800.0 + Mp, 2) - pow(800.0, 2);
    double xF = 0.2;
    while (infile >> var[0] >> var[1] >> value >> error[0]){
      Variable[Npt][0] = var[0];
      Variable[Npt][1] = var[1];
      Variable[Npt][2] = asinh(sqrt(var[3]) / var[0] * xF / 2.0);
      Variable[Npt][3] = var[3];
      Value[Npt] = value * 1.0e-10 / pow(0.197327, 2);
      Errors[Npt][0] = error[0] * 1.0e-10 / pow(0.197327, 2);
      Errors[Npt][1] = 0.10 * Value[Npt];
      if (var[0] < 6.0)
	Errors[Npt][1] = 0.20 * Value[Npt];
      Observable[Npt] = obs;
      Npt++;
    }
  }  
  infile.close();
  return 0;
}




#endif
