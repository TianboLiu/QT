#include "theory.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TString.h"

double Q[200], QT[200], ds[200], Estat[200], Esyst[200];
double s = pow(38.8, 2);
double xF = 0.1;
double dPA = 1.0;//proton
double dPB = 29.0 / 64.0;//Cu
int Npoints = 0;

int LoadData(const double Qvalue, const double QTmax, const char * file){
  ifstream infile(file);
  char ltmp[300];
  for (int i = 0; i < 15; i++)//skiprows
    infile.getline(ltmp, 300);
  int Npt = 0;
  while (infile >> Q[Npt] >> QT[Npt] >> ds[Npt] >> Estat[Npt]){
    if (abs(Q[Npt] - Qvalue) > 0.1) continue;
    if (QT[Npt] > QTmax) continue;
    //Esyst[Npt] = ds[Npt] * 0.15;
    Esyst[Npt] = 0;
    Npt++;
  }
  infile.close();
  return Npt;
}

double Chi2(const double * par){
  for (int i = 0; i < 3; i++)
    Parameters[i] = par[i];
  double theory;
  double sum = 0;
  double var[6] = {0.0, 0.0, 0.0, s, dPA, dPB};
  double unit = pow(1.0e13 / 0.197327, 2);
  for (int i = 0; i < Npoints; i++){
    var[0] = Q[i];
    var[1] = QT[i];
    var[2] = asinh(sqrt(s) / Q[i] * xF / 2.0);
    theory = dsigma_DY(var) * pow(Q[i], 2) * 2.0 * log((Q[i] + 0.5) / (Q[i] - 0.5)) / M_PI;
    if (Q[i] > 11.5 && Q[i] < 13.5)
      theory = theory / log((Q[i] + 0.5) / (Q[i] - 0.5)) * log((Q[i] + 1.0) / (Q[i] - 1.0));
    if (Q[i] > 13.5)
      theory = theory / log((Q[i] + 0.5) / (Q[i] - 0.5)) * log((Q[i] + 2.25) / (Q[i] - 2.25));
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
    cout << "./fit_E605 <Q>" << endl;
    return 0;
  }

  TString filename = "path/DY/DY.E605.list";

  F_TMD = & F_TMD_simple;
  double init[3] = {1.0, 2.0, 1.0};

  double Qvalue = atof(argv[1]);
  double Tfrac = 0.0;
  double QTmax = 0.0;
  double central[3], cov[9];
  double chi2;
  int npt = 5;
 

  FILE * fs = fopen("fs.dat", "w");
  fprintf(fs, "E605\n");
  fprintf(fs, "s:\t%.4f\n", s);
  fprintf(fs, "xF:\t%.4f\n", xF);
  fprintf(fs, "Q:\t%.4f\n\n", Qvalue);
  fprintf(fs, "QT/Q\tQT\tNpt\tChi2\tN[0]\tp[1]\tc[2]\tcov00\tcov01\tcov02\tcov10\tcov11\tcov12\tcov20\tcov21\tcov22\n");

  int lf = 0;
  double Tflist[3] = {0.25, 0.5, 1.0};
  while (lf < 3){
    Tfrac = Tflist[lf++];
    printf("%.2f\r", Tfrac);
    QTmax = Tfrac * Qvalue;
    Npoints = LoadData(Qvalue, QTmax, filename.Data());
    //if (! (Npoints > npt)) continue;
    npt = Npoints;
    chi2 = Minimize(init, central, cov) / (Npoints - 3);
    cout << chi2 << endl;
    fprintf(fs, "%.2f\t%.2E\t%d\t%.2E\t%.2E\t%.2E\t%.2E\t%.2E\t%.2E\t%.2E\t%.2E\t%.2E\t%.2E\t%.2E\t%.2E\t%.2E\n",
	    Tfrac, QTmax, Npoints, chi2,
	    central[0], central[1], central[2],
	    cov[0], cov[1], cov[2],
	    cov[3], cov[4], cov[5],
	    cov[6], cov[7], cov[8]);
  }
  fclose(fs);
  
  return 0;
}
