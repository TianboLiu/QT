#include "tmdevol3.h"
#include "load.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLIntegrator.h"

using namespace std;

int NPar = 5;
double Parameters[5] = {0.1, 0.1, 0.1, 0.1, 0.1};
double ParametersErrors[5] = {0.1, 0.1, 0.1, 0.1, 0.1};
double dProtonA = 1.0;
double dProtonB = 1.0;

int GetTMDs(double * tmds, const double dProton, const double x, const double b_T, const double Q){
  tmds[2] = dProton * TMDEVOL::F_TMD(2, x, b_T, Q) + (1.0 - dProton) * TMDEVOL::F_TMD(1, x, b_T, Q);
  tmds[1] = dProton * TMDEVOL::F_TMD(1, x, b_T, Q) + (1.0 - dProton) * TMDEVOL::F_TMD(2, x, b_T, Q);
  tmds[2+6] = dProton * TMDEVOL::F_TMD(-2, x, b_T, Q) + (1.0 - dProton) * TMDEVOL::F_TMD(-1, x, b_T, Q);
  tmds[1+6] = dProton * TMDEVOL::F_TMD(-1, x, b_T, Q) + (1.0 - dProton) * TMDEVOL::F_TMD(-2, x, b_T, Q);
  tmds[0] = TMDEVOL::F_TMD(21, x, b_T, Q);
  for (int i = 3; i <= 6; i++){
    tmds[i] = TMDEVOL::F_TMD(i, x, b_T, Q);
    tmds[i+6] = TMDEVOL::F_TMD(-i, x, b_T, Q);
  }
  return 0;
}

int dF_TMD(const int flavor, const double x, const double b_T, const double Q){
  double step = 1.0e-3;
  double value = (TMDEVOL::F_TMD(flavor, x, b_T + step, Q) - TMDEVOL::F_TMD(flavor, x, b_T - step, Q)) / (2.0 * step);
  return value;
}

double F_NP_1(const int flavor, const double x, const double b_T){
  return 1.0;
}

double F_NP_N(const int flavor, const double x, const double b_T){
  return Parameters[0];
}

const double alpha_EM_0 = 1.0 / 137.0;
const double e_u = 2.0 / 3.0;
const double e_d = -1.0 / 3.0;
double dsigma_DY_integrand(const double ATan_b_T, void * par){
  double * param = (double *) par;
  double Q = param[0];
  double QT = param[1];
  double y = param[2];
  double s = param[3];
  double xA = Q / sqrt(s) * exp(y);
  double xB = Q / sqrt(s) * exp(-y);
  double b_T = tan(ATan_b_T);
  double sigma0 = 4.0 * pow(M_PI, 2) * pow(alpha_EM_0, 2) / (9.0 * s * Q * Q);
  double tmdA[13], tmdB[13];
  GetTMDs(tmdA, dProtonA, xA, b_T, Q);
  GetTMDs(tmdB, dProtonB, xB, b_T, Q);
  double convol =
    pow(e_u, 2) * (tmdA[2] * tmdB[2+6] + tmdA[2+6] * tmdB[2]
		   + tmdA[4] * tmdB[4+6] + tmdA[4+6] * tmdB[4])
    +
    pow(e_d, 2) * (tmdA[1] * tmdB[1+6] + tmdA[1+6] * tmdB[1]
		   + tmdA[3] * tmdB[3+6] + tmdA[3+6] * tmdB[3]
		   + tmdA[5] * tmdB[5+6] + tmdA[5+6] * tmdB[5]);
  double result = b_T * ROOT::Math::cyl_bessel_j(0, b_T * QT) / (2.0 * M_PI) * sigma0 * convol;
  return result / pow(cos(ATan_b_T), 2);
}

double dsigma_DY(const double Q, const double QT, const double y, const double s){
  double par[4] = {Q, QT, y, s};
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-6);
  ig.SetFunction(&dsigma_DY_integrand, par);
  double result = ig.Integral(1e-6, M_PI/2.0 - 1e-2);
  //double result = ig.Integral(0.1, 0.2);
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
    if (Observable[i] == 0){//E288
      dProtonA = 1.0;
      dProtonB = 78.0 / 195.0;
      theory = dsigma_DY(Variable[i][0], Variable[i][1], Variable[i][2], Variable[i][3])
	* pow(Variable[i][0], 2) * 2.0 * log((Variable[i][0] + 0.5) / (Variable[i][0] - 0.5)) / M_PI;
      sum += pow(theory - Value[i], 2) / (pow(Errors[i][0], 2) + pow(Errors[i][1], 2));
    }
    else if (Observable[i] == 1){
      dProtonA = 1.0;
      dProtonB = 29.0 / 64.0;
      double dQ = 0.5;
      if (Variable[i][0] > 11.5) dQ = 1.0;
      if (Variable[i][0] > 13.5) dQ = (18.0 - 13.5) / 2.0;
      theory = dsigma_DY(Variable[i][0], Variable[i][1], Variable[i][2], Variable[i][3])
	* pow(Variable[i][0], 2) * 2.0 * log((Variable[i][0] + dQ) / (Variable[i][0] - dQ)) / M_PI;
      sum += pow(theory - Value[i], 2) / (pow(Errors[i][0], 2) + pow(Errors[i][1], 2));
    }
    else if (Observable[i] == 2){
      dProtonA = 1.0;
      dProtonB = 1.0 / 2.0;
      theory = dsigma_DY(Variable[i][0], Variable[i][1], Variable[i][2], Variable[i][3])
	* pow(Variable[i][0], 2) * 2.0 * log((Variable[i][0] + 0.5) / (Variable[i][0] - 0.5)) / M_PI;
      sum += pow(theory - Value[i], 2) / (pow(Errors[i][0], 2) + pow(Errors[i][1], 2));
    }    
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
    min->SetVariable(i, "p", init[i], 1.0e-4);
  }
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
    return 0;
  }


  int task = atoi(argv[1]);

  if (task == 0){
    LoadData_DY("path/Data/DY/DY.E288_200.list", "E288_200");
    LoadData_DY("path/Data/DY/DY.E288_300.list", "E288_300");
    LoadData_DY("path/Data/DY/DY.E288_400.list", "E288_400");
    LoadData_DY("path/Data/DY/DY.E605.list", "E605");
    LoadData_DY("path/Data/DY/DY.E772.list", "E772");

    double Q = atof(argv[2]);
    SetFlagQ(Q, 1e-3);
    SetFlagQT(2.0);

    int npt = CountPoints();
    cout << "Npoints: " << npt << " " << Npt << endl;

    TMDEVOL::Initialize();
    TMDEVOL::F_NP = & F_NP_1;
    NPar = 0;

    //for (int i = 0; i < Npt; i++){
    //  if (FlagQ[i])
    //	cout << Variable[i][0] << " " << Variable[i][1] << " " << Variable[i][2] << " " << Variable[i][3] << endl;
    //}

    cout << Chi2(Parameters) / npt << endl;
    
  }





  
  return 0;
}

