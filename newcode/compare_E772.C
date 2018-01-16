#include "theory.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TString.h"
#include "TMatrixDEigen.h"

double Q[200], QT[200], ds[200], Estat[200], Esyst[200];
double s = pow(38.8, 2);
double xF = 0.2;
double dPA = 1.0;//proton
double dPB = 1.0 / 2.0;//d
int Npoints = 0;
double unit = 1e-10 / pow(0.197327, 2);
double pars[7][3];
double Qvalue = 0;

int LoadData(const double Qvalue, const double QTmax, const char * file){
  FILE * fs = fopen("data_in.dat", "w");
  ifstream infile(file);
  char ltmp[300];
  int flag = 0;
  for (int i = 0; i < 15; i++)//skiprows
    infile.getline(ltmp, 300);
  int Npt = 0;
  while (infile >> Q[Npt] >> QT[Npt] >> ds[Npt] >> Estat[Npt]){
    if (abs(Q[Npt] - Qvalue) > 0.1) continue;
    Esyst[Npt] = ds[Npt] * 0.15;
    if (QT[Npt] < QTmax) {
      flag = 1;
      fprintf(fs, "%.4E\t%.4E\t%.4E\t%.4E\t%d\n",
	      QT[Npt], ds[Npt] * unit, Estat[Npt] * unit, Esyst[Npt] * unit, flag);
      Npt++;
    }
    else {
      flag = 0;
      fprintf(fs, "%.4E\t%.4E\t%.4E\t%.4E\t%d\n",
	      QT[Npt], ds[Npt] * unit, Estat[Npt] * unit, Esyst[Npt] * unit, flag);
    }
  }
  infile.close();
  fclose(fs);
  return Npt;
}

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

int Curve(){
  FILE * fc = fopen("fc.dat", "w");
  double theory[7];
  double error;
  double var[6] = {Qvalue, 0.0, asinh(sqrt(s) / Qvalue * xF / 2.0), s, dPA, dPB};
  for (int ix = 0; ix < 50; ix++){
    var[1] = 0.05 + 0.1 * ix;
    for (int i = 0; i < 7; i++){
      for (int j = 0; j < 3; j++)
	Parameters[j] = pars[i][j];
      theory[i] = dsigma_DY(var) * pow(Qvalue, 2) * 2.0 * log((Qvalue + 0.5) / (Qvalue - 0.5)) / M_PI;
    }
    error = sqrt(pow(theory[1] - theory[2], 2) + pow(theory[3] - theory[4], 2) + pow(theory[5] - theory[6], 2)) / 2.0;
    fprintf(fc, "%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\n",
	    var[1], theory[0], error, theory[1], theory[2], theory[3], theory[4], theory[5], theory[6]);
  }
  return 0;
}


int main(const int argc, const char * argv[]){

  if (argc < 4){
    cout << "./compare_E772 <Q> <QT> <fitfile>" << endl;
    return 0;
  }

  TString filename = "path/DY/DY.E772.list";

  F_TMD = & F_TMD_simple;

  Qvalue = atof(argv[1]);
  double QTmax = atof(argv[2]);

  LoadData(Qvalue, QTmax, filename.Data());
  
  GetParameters(argv[3]);
  Curve();
  
  return 0;
}
