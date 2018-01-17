#include "theory.h"
#include "tmdevol3.h"

#include "TMatrixDEigen.h"

double pars[7][3];

int GetParameters(const char * file){
  ifstream fp(file);
  char tmp[300];
  for (int i = 0; i < 8; i++)
    fp.getline(tmp, 300);
  fp >> pars[0][0] >> pars[0][1] >> pars[0][2];
  fp.getline(tmp, 300);
  TMatrixD cov(3,3);
  fp >> cov(0,0) >> cov(0,1) >> cov(0,2);
  fp >> cov(1,0) >> cov(1,1) >> cov(1,2);
  fp >> cov(2,0) >> cov(2,1) >> cov(2,2);
  TMatrixDEigen eigen(cov);
  TMatrixD val = eigen.GetEigenValues();
  TMatrixD vec = eigen.GetEigenVectors();
  pars[1][0] = pars[0][0] + val(0,0) * vec(0,0);
  pars[1][1] = pars[0][1] + val(0,0) * vec(1,0);
  pars[1][2] = pars[0][2] + val(0,0) * vec(2,0);
  pars[2][0] = pars[0][0] - val(0,0) * vec(0,0);
  pars[2][1] = pars[0][1] - val(0,0) * vec(1,0);
  pars[2][2] = pars[0][2] - val(0,0) * vec(2,0);
  pars[3][0] = pars[0][0] + val(1,1) * vec(0,1);
  pars[3][1] = pars[0][1] + val(1,1) * vec(1,1);
  pars[3][2] = pars[0][2] + val(1,1) * vec(2,1);
  pars[4][0] = pars[0][0] - val(1,1) * vec(0,1);
  pars[4][1] = pars[0][1] - val(1,1) * vec(1,1);
  pars[4][2] = pars[0][2] - val(1,1) * vec(2,1);
  pars[5][0] = pars[0][0] + val(2,2) * vec(0,2);
  pars[5][1] = pars[0][1] + val(2,2) * vec(1,2);
  pars[5][2] = pars[0][2] + val(2,2) * vec(2,2);
  pars[6][0] = pars[0][0] - val(2,2) * vec(0,2);
  pars[6][1] = pars[0][1] - val(2,2) * vec(1,2);
  pars[6][2] = pars[0][2] - val(2,2) * vec(2,2);
  return 0;
}

double F_TMD_evol(const int flavor, const double x, const double b_T, const double Q0, const double Q){
  return F_TMD_simple(flavor, x, b_T, Q0) * exp(TMDEVOL::S(b_T, Q0) - TMDEVOL::S(b_T, Q));
}

int InputCurve(const double x, const double Q0){
  FILE * fQ0 = fopen("fQ0.dat", "w");
  double fu[7];
  double error;
  double b;
  for (int ib = 0; ib < 50; ib++){
    b = 0.05 + 0.1 * ib;
    for (int i = 0; i < 7; i++){
      for (int j = 0; j < 3; j++)
	Parameters[j] = pars[i][j];
      fu[i] = F_TMD_simple(2, x, b, Q0);
    }
    error = sqrt(pow(fu[1] - fu[2], 2) + pow(fu[3] - fu[4], 2) + pow(fu[5] - fu[6], 2)) / 2.0;
    fprintf(fQ0, "%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\n",
	    b, fu[0], error, fu[1], fu[2], fu[3], fu[4], fu[5], fu[6]);
  }
  fclose(fQ0);
  return 0;
}

int EvolvedCurve(const double x, const double b, const double Q0){
  FILE * fQ = fopen("fQ.dat", "w");
  double fu[7];
  double error;
  double Q;
  for (int iQ = 0; iQ < 87; iQ++){
    Q = 4.5 + 1.0 * iQ;
    cout << Q << endl;
    for (int i = 0; i < 7; i++){
      for (int j = 0; j < 3; j++)
	Parameters[j] = pars[i][j];
      fu[i] = F_TMD_evol(2, x, b, Q0, Q);
    }
    error = sqrt(pow(fu[1] - fu[2], 2) + pow(fu[3] - fu[4], 2) + pow(fu[5] - fu[6], 2)) / 2.0;
    fprintf(fQ, "%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\n",
	    Q, fu[0], error, fu[1], fu[2], fu[3], fu[4], fu[5], fu[6]);
  }
  fclose(fQ);
  return 0;
}

int main(const int argc, const char * argv[]){

  if (argc < 2){
    cout << "./calculate <opt>" << endl;
    return 0;
  }

  if (strcmp(argv[1], "fitcurve") == 0){
    if (argc < 5){
      cout << "./calculate fitcurve <x> <Q0> <fitfile>" << endl;
      return 0;
    }
    double x = atof(argv[2]);
    double Q0 = atof(argv[3]);
    GetParameters(argv[4]);
    InputCurve(x, Q0);
  }

  if (strcmp(argv[1], "evolcurve") == 0){
    if (argc < 6){
      cout << "./calculate fitcurve <x> <b> <Q0> <fitfile>" << endl;
      return 0;
    }
    TMDEVOL::Initialize();
    double x = atof(argv[2]);
    double b = atof(argv[3]);
    double Q0 = atof(argv[4]);
    GetParameters(argv[5]);
    EvolvedCurve(x, b, Q0);
  }
  
  return 0;
}
