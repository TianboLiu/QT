#include "theory.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TString.h"
#include "TMatrixDEigen.h"

double QT[200], ds[200], Estat[200], Esyst[200];
double s = pow(1960.0, 2);
double dPA = 1.0;//proton
double dPB = -1.0;//antiproton
double unit = 255.8 * 1e-10 / pow(0.197327, 2);
double pars[7][3];
int Npoints = 0;

int LoadData(const double QTmax, const char * file){
  FILE * fs = fopen("data_in.dat", "w");
  ifstream infile(file);
  char ltmp[300];
  int flag = 0;
  for (int i = 0; i < 20; i++)//skiprows
    infile.getline(ltmp, 300);
  int Npt = 0;
  while (infile >> QT[Npt] >> ds[Npt] >> Estat[Npt] >> Esyst[Npt]){
    if (QT[Npt] < QTmax){
      flag = 1;
      fprintf(fs, "%.4E\t%.4E\t%.4E\t%.4E\t%d\n",
	      QT[Npt], ds[Npt] * unit, Estat[Npt] * unit, 0.0, flag);
      Npt++;
    }
    else {
      flag = 0;
      fprintf(fs, "%.4E\t%.4E\t%.4E\t%.4E\t%d\n",
	      QT[Npt], ds[Npt] * unit, Estat[Npt] * unit, 0.0, flag);
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

double Integrand(const double * var, const double * par){
  double Q = var[0];
  double y = var[1];
  double QT = par[0];
  double para[6] = {Q, QT, y, s, dPA, dPB};
  double xs = dsigma_Z(para) + dsigma_gZ(para);
  return 2.0 * Q * 2.0 * QT * xs;
}

double Theory(const double QT){
  double xl[2] = {66.0, -3.39};
  double xu[2] = {116.0, 3.39};
  ROOT::Math::WrappedParamFunction<> wf(&Integrand, 2, 1);
  wf.SetParameters(&QT);
  ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE, 0.0, 1e-5, 50000);
  ig.SetFunction(wf);
  double result = ig.Integral(xl, xu);
  return result;
}
  
int Curve(){
  FILE * fc = fopen("fc.dat", "w");
  double theory[7];
  double error;
  double QT;
  for (int ix = 0; ix < 50; ix++){
    QT = 0.05 + 0.5 * ix;
    for (int i = 0; i < 7; i++){
      for (int j = 0; j < 3; j++)
	Parameters[j] = pars[i][j];
      theory[i] = Theory(QT);
    }
    error = sqrt(pow(theory[1] - theory[2], 2) + pow(theory[3] - theory[4], 2) + pow(theory[5] - theory[6], 2)) / 2.0;
    fprintf(fc, "%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\n",
	    QT, theory[0], error, theory[1], theory[2], theory[3], theory[4], theory[5], theory[6]);
  }
  fclose(fc);
  return 0;
}

int main(const int argc, const char * argv[]){

  if (argc < 3){
    cout << "./fit_D0II <QT> <fitfile>" << endl;
    return 0;
  }

  TString filename = "path/DY/DY.D0_RunII.list";

  F_TMD = & F_TMD_simple;

  double QTmax = atof(argv[1]);

  LoadData(QTmax, filename.Data());
    
  GetParameters(argv[2]);
  Curve();
  
  return 0;
}
