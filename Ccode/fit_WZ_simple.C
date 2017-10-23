#include "LHAPDF/LHAPDF.h"

#include "tmdevol3.h"
#include "load.h"

#include "TString.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLIntegrator.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraph.h"
#include "TH1D.h"

using namespace std;
LHAPDF::PDF * xpdf = LHAPDF::mkPDF("CJ15lo", 0);

int NPar = 5;
double Parameters[5] = {1.0, 6.0 / 3.0, 1.0, 0.1, 0.1};
double ParametersErrors[5] = {0.0, 0.1, 0.1, 0.1, 0.1};
double dProtonA = 1.0;
double dProtonB = 1.0;

double Q0;

double (* F_TMD)(const int flavor, const double x, const double b_T, const double Q);

double F_TMD_simple(const int flavor, const double x, const double b_T, const double Q){
  double col = xpdf->xfxQ(flavor, x, Q) / x;
  double trans = exp(-Parameters[2] * pow(b_T, Parameters[1]));
  return Parameters[0] * col * trans;
}

double F_TMD_evol(const int flavor, const double x, const double b_T, const double Q){
  return F_TMD_simple(flavor, x, b_T, Q0) * exp(TMDEVOL::S(b_T, Q0) - TMDEVOL::S(b_T, Q));
}

int Exchange(double &a, double &b){
  double c = a;
  a = b;
  b = c;
  return 0;
}

int GetTMDs(double * tmds, const double dProton, const double x, const double b_T, const double Q){
  tmds[2] = abs(dProton) * F_TMD(2, x, b_T, Q) + (1.0 - abs(dProton)) * F_TMD(1, x, b_T, Q);
  tmds[1] = abs(dProton) * F_TMD(1, x, b_T, Q) + (1.0 - abs(dProton)) * F_TMD(2, x, b_T, Q);
  tmds[2+6] = abs(dProton) * F_TMD(-2, x, b_T, Q) + (1.0 - abs(dProton)) * F_TMD(-1, x, b_T, Q);
  tmds[1+6] = abs(dProton) * F_TMD(-1, x, b_T, Q) + (1.0 - abs(dProton)) * F_TMD(-2, x, b_T, Q);
  tmds[0] = F_TMD(21, x, b_T, Q);
  for (int i = 3; i <= 6; i++){
    tmds[i] = F_TMD(i, x, b_T, Q);
    tmds[i+6] = F_TMD(-i, x, b_T, Q);
  }
  if (dProton < 0){
    for (int i = 1; i <= 6; i++){
      Exchange(tmds[i], tmds[i+6]);
    }
  }
  return 0;
}


const double alpha_EM_0 = 1.0 / 137.0;
const double alpha_EM_Z = 1.0 / 128/0;
const double e_u = 2.0 / 3.0;
const double e_d = -1.0 / 3.0;
const double sinW2 = 0.2313;
const double V_u = 0.5 - 2.0 * e_u * sinW2;
const double V_d = -0.5 - 2.0 * e_d * sinW2;
const double A_u = 0.5;
const double A_d = -0.5;

double dsigma_Z_integrand(const double b_T, void * par){
  double * param = (double *) par;
  double Q = param[0];
  double QT = param[1];
  double y = param[2];
  double s = param[3];
  double xA = Q / sqrt(s) * exp(y);
  double xB = Q / sqrt(s) * exp(-y);
  if (xA >= 1.0 || xB >= 1.0) return 0;
  double sigma0 = pow(M_PI, 2) * alpha_EM_Z / (3.0 * s * sinW2 * (1.0 - sinW2));
  double pdfA[13], pdfB[13];
  GetTMDs(pdfA, dProtonA, xA, b_T, Q);
  GetTMDs(pdfB, dProtonB, xB, b_T, Q);
  double convol =
    (pow(V_u, 2) + pow(A_u, 2)) * (pdfA[2] * pdfB[2+6] + pdfA[2+6] * pdfB[2]
		   + pdfA[4] * pdfB[4+6] + pdfA[4+6] * pdfB[4])
    +
    (pow(V_d, 2) + pow(A_d, 2)) * (pdfA[1] * pdfB[1+6] + pdfA[1+6] * pdfB[1]
		   + pdfA[3] * pdfB[3+6] + pdfA[3+6] * pdfB[3]
		   + pdfA[5] * pdfB[5+6] + pdfA[5+6] * pdfB[5]);
  double result = b_T * ROOT::Math::cyl_bessel_j(0, b_T * QT) / (2.0 * M_PI) * sigma0 * convol * 0.03366;
  return result;
}

double dsigma_Z(const double Q, const double QT, const double y, const double s){
  double par[4] = {Q, QT, y, s};
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-6);
  ig.SetFunction(&dsigma_Z_integrand, par);
  double result = ig.Integral(1e-6, 70.0);
  return result;
}

double Z_integrand_y(const double y, void * par){
  double * param = (double *) par;
  double result = dsigma_Z(param[0], param[1], y, param[3]);
  return result;
}

double dsigma_Z_QT(const double Q, const double QT, const double s){
  double ylim = log(sqrt(s) / Q);
  double par[4] = {Q, QT, 0.0, s};
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-6);
  ig.SetFunction(&Z_integrand_y, par);
  double result = ig.Integral(-ylim, ylim);
  return result;
}

double Chi2(const double * par = Parameters){
  for (int i = 0; i < NPar; i++){
    Parameters[i] = par[i];
  }
  double theory;
  double sum = 0.0;
  for (int i = 0; i < Npt; i++){
    if (!(FlagQ[i] && FlagQT[i])) continue;
    if (Observable[i] == 10){//p pbar -> Z
      dProtonA = 1.0;
      dProtonB = -1.0;
      theory = dsigma_Z_QT(Variable[i][0], Variable[i][1], Variable[i][3])
        * 2.0 * Varible[i][1];
      sum += pow(theory - Value[i], 2) / (pow(Errors[i][0], 2) + pow(Errors[i][1], 2));
    }
  } 
  return sum;
}

int option = 0;
double Minimize(const int NPAR, const double * init){
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(100000);
  //min->SetPrecision(1.0e-14);
  min->SetTolerance(1.0e-3);
  min->SetPrintLevel(0);
  ROOT::Math::Functor f(&Chi2, NPAR);
  min->SetFunction(f);
  //for (int i = 0; i < NPAR; i++){
  //  min->SetVariable(i, "p", init[i], 1.0e-4);
  //}
  min->SetVariable(0, "N", init[0], 1e-4);
  if (option == 0) min->SetFixedVariable(0, "N", 1.0);
  min->SetLimitedVariable(1, "alpha", init[1], 1e-4, 0.01, 10.0);
  min->SetLowerLimitedVariable(2, "c", init[2], 1e-4, 0.001);
  min->Minimize();
  const double * xs = min->X();
  const double * es = min->Errors();
  for (int i = 0; i < NPAR; i++){
    Parameters[i] = xs[i];
    ParametersErrors[i] = es[i];
  }
  const double chi2 = min->MinValue();
  return chi2;
}

int CountPoints(){
  int n = 0;
  for (int i = 0; i < Npt; i++){
    if (FlagQ[i] && FlagQT[i]) n++;
  }
  return n;
} 


int main(const int argc, const char * argv[]){

  if (argc < 2){
    cout << "./compare <dataset> <task>" <<  endl;
    return 0;
  }

  int dataset = atoi(argv[1]);
  int task = atoi(argv[2]);
  option = 1;

  if (data == 0){
    LoadData_Z("path/Data/DY/DY.CDF_RunI.list", "CDF_RunI");
    LoadData_Z("path/Data/DY/DY.CDF_RunII.list", "CDF_RunII");
    LoadData_Z("path/Data/DY/DY.D0_RunI.list", "D0_RunI");
    LoadData_Z("path/Data/DY/DY.D0_RunII.list", "D0_RunII");
  }
  else if (data == 1){
    LoadData_Z("path/Data/DY/DY.CDF_RunI.list", "CDF_RunI");
  }

  if (task == 0){//QT scan
    F_TMD = & F_TMD_simple;

    double Q = atof(argv[3]);
    SetFlagQ(Q, 1.0);
    double QT = 0.0;

    int npt = 0;
    printf("Q\tQT\tnpt\tchi2\tN\talpha\tc\tdN\tdalpha\tdc\n");
    while (QT < Q){
      QT += 0.1;
      SetFlagQT(QT);
      if (npt == CountPoints()) continue;
      npt = CountPoints();

      if (npt < 4) continue;

      NPar = 3;
      Parameters[0] = 1.0;
      Parameters[1] = 2.0;
      Parameters[2] = 1.0;
      double chi2 = 0.0;

      chi2 = Minimize(NPar, Parameters) / (npt - 3);

      printf("%.1f\t%.1f\t%d\t%.2f\t%.2E\t%.2E\t%.2E\t%.1E\t%.1E\t%.1E\n",
	     Q, QT, npt, chi2, Parameters[0], Parameters[1], Parameters[2], ParametersErrors[0], ParametersErrors[1], ParametersErrors[2]);
    }
  }


  return 0;
}

