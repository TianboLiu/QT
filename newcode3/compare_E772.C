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

int LoadData(const double Qvalue, const char * file){
  FILE * fs = fopen("data.dat", "w");
  ifstream infile(file);
  char ltmp[300];
  for (int i = 0; i < 15; i++)//skiprows
    infile.getline(ltmp, 300);
  int Npt = 0;
  while (infile >> Q[Npt] >> QT[Npt] >> ds[Npt] >> Estat[Npt]){
    if (abs(Q[Npt] - Qvalue) > 0.1) continue;
    //Esyst[Npt] = ds[Npt] * 0.15;
    Esyst[Npt] = 0;
    fprintf(fs, "%.4E\t%.4E\t%.4E\t%.4E\n",
	      QT[Npt], ds[Npt] * unit, Estat[Npt] * unit, Esyst[Npt] * unit);
    Npt++;
  }
  infile.close();
  fclose(fs);
  return Npt;
}

int GetParameters(const char * file, const int line){
  ifstream fp(file);
  char tmp[300];
  for (int i = 0; i < 6 + line; i++)
    fp.getline(tmp, 300);
  double temp;
  TMatrixD cov(3,3);
  fp >> temp >> temp >> temp >> temp
     >> pars[0][0] >> pars[0][1] >> pars[0][2]
     >> cov(0,0) >> cov(0,1) >> cov(0,2)
     >> cov(1,0) >> cov(1,1) >> cov(1,2)
     >> cov(2,0) >> cov(2,1) >> cov(2,2);
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
  for (int i = 0; i < 7; i++){
    while (pars[i][0] < 0 || pars[i][1] < 0 || pars[i][2] < 0){
      pars[i][0] = 0.1 * pars[0][0] + 0.9 * pars[i][0];
      pars[i][1] = 0.1 * pars[0][1] + 0.9 * pars[i][1];
      pars[i][2] = 0.1 * pars[0][2] + 0.9 * pars[i][2];
    }
  }
  return 0;
}

int Curve(const char * filename, const double Tf){
  FILE * fc = fopen(filename, "w");
  double theory[7];
  double error;
  double var[6] = {Qvalue, 0.0, asinh(sqrt(s) / Qvalue * xF / 2.0), s, dPA, dPB};
  for (int ix = 0; ix < 50; ix++){
    var[1] = 0.05 + 0.1 * ix;
    if (var[1] > Tf * Qvalue) break;
    for (int i = 0; i < 7; i++){
      for (int j = 0; j < 3; j++)
	Parameters[j] = pars[i][j];
      theory[i] = dsigma_DY(var) * pow(Qvalue, 2) * 2.0 * log((Qvalue + 0.5) / (Qvalue - 0.5)) / M_PI;
    }
    error = sqrt(pow(theory[1] - theory[2], 2) + pow(theory[3] - theory[4], 2) + pow(theory[5] - theory[6], 2)) / 2.0;
    fprintf(fc, "%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\n",
	    var[1], theory[0], error, theory[1], theory[2], theory[3], theory[4], theory[5], theory[6]);
  }
  fclose(fc);
  return 0;
}


int main(const int argc, const char * argv[]){

  if (argc < 3){
    cout << "./compare_E772 <Q> <fitfile>" << endl;
    return 0;
  }

  TString filename = "path/DY/DY.E772.list";

  F_TMD = & F_TMD_simple;

  Qvalue = atof(argv[1]);

  TString curve1 = "fc0.dat";
  TString curve2 = "fc1.dat";
  TString curve3 = "fc2.dat";
  LoadData(Qvalue, filename.Data());

  GetParameters(argv[2], 0);
  Curve(curve1.Data(), 0.25);

  GetParameters(argv[2], 1);
  Curve(curve2.Data(), 0.5);
  
  GetParameters(argv[2], 2);
  Curve(curve3.Data(), 1.0);
  
  return 0;
}
