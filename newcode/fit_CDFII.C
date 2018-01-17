#include "theory.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TString.h"

double QT[200], ds[200], Estat[200], Esyst[200];
double s = pow(1960.0, 2);
double dPA = 1.0;//proton
double dPB = -1.0;//antiproton
int Npoints = 0;

int LoadData(const double QTmax, const char * file){
  ifstream infile(file);
  char ltmp[300];
  for (int i = 0; i < 20; i++)//skiprows
    infile.getline(ltmp, 300);
  int Npt = 0;
  while (infile >> QT[Npt] >> ds[Npt] >> Estat[Npt] >> Esyst[Npt]){
    if (QT[Npt] > QTmax) continue;
    Npt++;
  }
  infile.close();
  return Npt;
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
  double xl[2] = {66.0, -4.2};
  double xu[2] = {116.0, 4.2};
  ROOT::Math::WrappedParamFunction<> wf(&Integrand, 2, 1);
  wf.SetParameters(&QT);
  ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE, 0.0, 1e-4, 10000);
  ig.SetFunction(wf);
  double result = ig.Integral(xl, xu);
  return result;
}
  
double Chi2(const double * par){
  for (int i = 0; i < 3; i++)
    Parameters[i] = par[i];
  double theory;
  double sum = 0;
  double unit = 1e-10 / pow(0.197327, 2);
  for (int i = 0; i < Npoints; i++){
    theory = Theory(QT[i]);
    sum += pow(theory - ds[i] * unit, 2) / (pow(Estat[i] * unit, 2) + pow(Esyst[i] * unit, 2));
  } 
  return sum;
}

double Minimize(const double * init, double * central, double * cov){
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(1.0e-3);
  min->SetPrintLevel(1);
  ROOT::Math::Functor f(&Chi2, 3);
  min->SetFunction(f);
  min->SetLowerLimitedVariable(0, "N", init[0], 1e-4, 0.0);
  min->SetLimitedVariable(1, "a", init[1], 1e-4, 0.0, 5.0);
  min->SetLimitedVariable(2, "c", init[2], 1e-4, 0.0, 5.0);
  min->Minimize();
  const double * xs = min->X();
  for (int i = 0; i < 3; i++)
    central[i] = xs[i];
  min->GetCovMatrix(cov);
  return min->MinValue();
}

int main(const int argc, const char * argv[]){

  if (argc < 2){
    cout << "./fit_CDFII <QT>" << endl;
    return 0;
  }

  TString filename = "path/DY/DY.CDF_RunII.list";

  F_TMD = & F_TMD_simple;
  double init[3] = {1.0, 2.0, 1.0};

  double QTmax = atof(argv[1]);

  Npoints = LoadData(QTmax, filename.Data());
  cout << Npoints << endl;
  
  double central[3], cov[9];
  double chi2 = Minimize(init, central, cov) / (Npoints - 3);

  FILE * fs = fopen("fs.dat", "w");
  fprintf(fs, "CDF RunII\n");
  fprintf(fs, "s:\t%.4f\n", s);
  fprintf(fs, "y:\t integrated\n");
  fprintf(fs, "Q:\t 66-116\n\n");
  fprintf(fs, "%d\t%.2f\n", Npoints, chi2);
  fprintf(fs, "QT:\t%.4f\n\n", QTmax);
  fprintf(fs, "%.2E\t%.2E\t%.2E\n\n", central[0], central[1], central[2]);
  fprintf(fs, "%.2E\t%.2E\t%.2E\n%.2E\t%.2E\t%.2E\n%.2E\t%.2E\t%.2E\n",
	  cov[0], cov[1], cov[2],
	  cov[3], cov[4], cov[5],
	  cov[6], cov[7], cov[8]);
  fclose(fs);
  
  return 0;
}
