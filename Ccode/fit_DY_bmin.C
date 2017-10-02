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
double Parameters[5] = {1.0, 1.0 / 3.0, 1.0, 0.1, 0.1};
double ParametersErrors[5] = {0.0, 0.1, 0.1, 0.1, 0.1};
double dProtonA = 1.0;
double dProtonB = 1.0;

double Q0;
double bmin, bmax;

double (* F_TMD)(const int flavor, const double x, const double b_T, const double Q);

double F_TMD_full(const int flavor, const double x, const double b_T, const double Q){
  double N = Parameters[0];
  if (b_T < bmin){
    return N * TMDEVOL::F_TMD(flavor, x, bmin, Q);
  }
  else if (b_T < bmax){
    return N * TMDEVOL::F_TMD(flavor, x, b_T, Q);
  }
  else{
    double fmax = TMDEVOL::F_TMD(flavor, x, bmax, Q);
    if (fmax == 0) return 0;
    double c = Parameters[2];
    //double c = - TMDEVOL::dF_TMD_$_db_T(flavor, x, bmax, Q) / fmax / (Parameters[1] * pow(bmax, Parameters[1] - 1.0));
    return N * fmax * exp(-c * (pow(b_T, Parameters[1]) - pow(bmax, Parameters[1])));
  }
}

double F_TMD_evol(const int flavor, const double x, const double b_T, const double Q){
  return F_TMD(flavor, x, b_T, Q0) * exp(TMDEVOL::S(b_T, Q0) - TMDEVOL::S(b_T, Q));
}
  
int GetPDFs(double * pdfs, const double dProton, const double x, const double Q){
  pdfs[2] = dProton * xpdf->xfxQ(2, x, Q) + (1.0 - dProton) * xpdf->xfxQ(1, x, Q);
  pdfs[1] = dProton * xpdf->xfxQ(1, x, Q) + (1.0 - dProton) * xpdf->xfxQ(2, x, Q);
  pdfs[2+6] = dProton * xpdf->xfxQ(-2, x, Q) + (1.0 - dProton) * xpdf->xfxQ(-1, x, Q);
  pdfs[1+6] = dProton * xpdf->xfxQ(-1, x, Q) + (1.0 - dProton) * xpdf->xfxQ(-2, x, Q);
  pdfs[0] = xpdf->xfxQ(21, x, Q);
  for (int i = 3; i <= 6; i++){
    pdfs[i] = xpdf->xfxQ(i, x, Q);
    pdfs[i+6] = xpdf->xfxQ(-i, x, Q);
  }
  for (int i = 0; i < 13; i++)
    pdfs[i] = pdfs[i] / x;
  return 0;  
}

int GetTMDs(double * tmds, const double dProton, const double x, const double b_T, const double Q){
  tmds[2] = dProton * F_TMD(2, x, b_T, Q) + (1.0 - dProton) * F_TMD(1, x, b_T, Q);
  tmds[1] = dProton * F_TMD(1, x, b_T, Q) + (1.0 - dProton) * F_TMD(2, x, b_T, Q);
  tmds[2+6] = dProton * F_TMD(-2, x, b_T, Q) + (1.0 - dProton) * F_TMD(-1, x, b_T, Q);
  tmds[1+6] = dProton * F_TMD(-1, x, b_T, Q) + (1.0 - dProton) * F_TMD(-2, x, b_T, Q);
  tmds[0] = F_TMD(21, x, b_T, Q);
  for (int i = 3; i <= 6; i++){
    tmds[i] = F_TMD(i, x, b_T, Q);
    tmds[i+6] = F_TMD(-i, x, b_T, Q);
  }
  return 0;
}

const double alpha_EM_0 = 1.0 / 137.0;
const double e_u = 2.0 / 3.0;
const double e_d = -1.0 / 3.0;
double dsigma_DY_integrand(const double b_T, void * par){
  double * param = (double *) par;
  double Q = param[0];
  double QT = param[1];
  double y = param[2];
  double s = param[3];
  double xA = Q / sqrt(s) * exp(y);
  double xB = Q / sqrt(s) * exp(-y);
  //double b_T = tan(ATan_b_T);
  double sigma0 = 4.0 * pow(M_PI, 2) * pow(alpha_EM_0, 2) / (9.0 * s * Q * Q);
  double pdfA[13], pdfB[13];
  GetTMDs(pdfA, dProtonA, xA, b_T, Q);
  GetTMDs(pdfB, dProtonB, xB, b_T, Q);
  double convol =
    pow(e_u, 2) * (pdfA[2] * pdfB[2+6] + pdfA[2+6] * pdfB[2]
		   + pdfA[4] * pdfB[4+6] + pdfA[4+6] * pdfB[4])
    +
    pow(e_d, 2) * (pdfA[1] * pdfB[1+6] + pdfA[1+6] * pdfB[1]
		   + pdfA[3] * pdfB[3+6] + pdfA[3+6] * pdfB[3]
		   + pdfA[5] * pdfB[5+6] + pdfA[5+6] * pdfB[5]);
  double result = b_T * ROOT::Math::cyl_bessel_j(0, b_T * QT) / (2.0 * M_PI) * sigma0 * convol;
  return result;// / pow(cos(ATan_b_T), 2);
}

double dsigma_DY(const double Q, const double QT, const double y, const double s){
  double par[4] = {Q, QT, y, s};
  //cout << Q << " " << QT << " " << y << " " << s << endl;
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-4);
  ig.SetFunction(&dsigma_DY_integrand, par);
  //double result = ig.Integral(1e-3, M_PI / 2.0 - 1e-6);
  double result = ig.Integral(1e-3, bmin) + ig.Integral(bmin, bmax) + ig.Integral(bmax, 100.0);
  //double result = ig.Integral(1e-3, 10.0);
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
  min->SetLowerLimitedVariable(1, "alpha", init[1], 1e-4, 1e-6);
  min->SetLowerLimitedVariable(2, "c", init[2], 1e-4, 1e-6);
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
  option = atoi(argv[2]);

  if (task == 0){
    TMDEVOL::Initialize();
    F_TMD = & F_TMD_full;
    bmax = 1.0;
    bmin = 0.1;

    //cout << F_TMD(2, 0.1, 2.5, 4.5) << endl;
    //cout << TMDEVOL::dF_TMD_$_db_T(2, 0.1, 1.0, 4.5) << endl;
    //cout << dsigma_DY(4.5, 0.3, 0.4, 19.4*19.4) << endl;
    //return 0;
    
    LoadData_DY("path/Data/DY/DY.E288_200.list", "E288_200");
    LoadData_DY("path/Data/DY/DY.E288_300.list", "E288_300");
    LoadData_DY("path/Data/DY/DY.E288_400.list", "E288_400");
    //LoadData_DY("path/Data/DY/DY.E605.list", "E605");
    //LoadData_DY("path/Data/DY/DY.E772.list", "E772");

    TString filename = "temp.dat";
    if (option == 0) filename = "results/fit_E288_bmin_0.dat";
    else if (option == 1) filename = "results/fit_E288_bmin_1.dat";
    else return 0;
    FILE * fs = fopen(filename.Data(), "w");
    
    fprintf(fs, "Q \t npt \t chi2 \t N \t alpha \t dN \t dalpha \n");

    double Q = 4.5;
    int npt = 0;
    double chi2 = 0;
    
    for (int j = 0; j < 10; j++){
      Q = 4.5 + j;
      Q0 = Q;
      bmin = 1.0 / Q;
      SetFlagQ(Q, 1e-3);
      SetFlagQT(112.0);
      npt = CountPoints();
      cout << "Npoints: " << npt << " " << Npt << endl;
      
      NPar = 3;
      Parameters[0] = 1.0;
      Parameters[1] = 0.5;
      Parameters[2] = 0.1;
      if (option == 0) chi2 = Minimize(NPar, Parameters) / (npt - 1);
      else if (option == 1) chi2 = Minimize(NPar, Parameters) / (npt - 2);
      else break;
 
      cout << Q << "  " << chi2 << "  " << bmin << ":  " << Parameters[0] << " " << Parameters[1] << Parameters[2] << endl;
      fprintf(fs, "%.1f \t %d \t %.2f \t %.3E \t %.3E \t %.3E \t %.3E \t %.3E \t %.3E\n",
	      Q, npt, chi2, Parameters[0], Parameters[1], Parameters[2], ParametersErrors[0], ParametersErrors[1], ParametersErrors[2]);
    }
    fclose(fs);
  }


  if (task == 10){
    F_TMD = &F_TMD_full;

    TMDEVOL::Initialize();
          
    double b1 = 0.1;
    double b2 = 0.2;
    double Q0;

    ifstream infile("results/fit_E288_bmin_0.dat");
    char tmp[300];
    infile.getline(tmp, 300);

    int colorlist[10] = {1, 4, 2, 6, 7, 1, 4, 2, 6, 7};
    TGraph g1(50);
    TGraph p1(1);
    TGraph g2(50);
    TGraph p2(50);
    double temp;

    TH1D * h0 = new TH1D("h0", "", 1, 3.5, 15.5);
    h0->SetStats(0);
    h0->SetMinimum(0);
    h0->SetMaximum(11.0);
    h0->GetXaxis()->SetTitle("Q (GeV)");
    h0->GetXaxis()->CenterTitle(true);
    h0->GetXaxis()->SetTitleSize(0.055);
    h0->GetXaxis()->SetTitleOffset(1.05);
    h0->GetXaxis()->SetLabelSize(0.055);
    h0->GetXaxis()->SetNdivisions(6, 5, 0);
    h0->GetYaxis()->SetTitle("f^u(x, b)");
    h0->GetYaxis()->CenterTitle(true);
    h0->GetYaxis()->SetTitleSize(0.055);
    h0->GetYaxis()->SetTitleOffset(1.05);
    h0->GetYaxis()->SetLabelSize(0.055);
    h0->GetYaxis()->SetNdivisions(6, 5, 0);

    TCanvas * c0 = new TCanvas("c0", "", 800, 1200);
    c0->Divide(1, 2);

    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);

    c0->cd(1);
    h0->DrawClone("AXIS");
    c0->cd(2);
    h0->DrawClone("AXIS");
    
    double Q[50];
    for (int i = 0; i < 50; i++)
      Q[i] = 4.5 + 0.2 * i;

    int j = 0;
    p1.SetMarkerStyle(20);
    p1.SetMarkerSize(0.5);
    p2.SetMarkerStyle(20);
    p2.SetMarkerSize(0.5);
    while(infile >> Q0 >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
      cout << Parameters[0] << " " << Parameters[1] << " " << Parameters[2] << endl;
      p1.SetPoint(0, Q0, F_TMD(2, 0.1, b1, Q0));
      p2.SetPoint(0, Q0, F_TMD(2, 0.1, b2, Q0));
      for (int i = 0; i < 50; i++){
	g1.SetPoint(i, Q[i], F_TMD(2, 0.1, b1, Q0) * exp(TMDEVOL::S(b1, Q0) - TMDEVOL::S(b1, Q[i])));
	g2.SetPoint(i, Q[i], F_TMD(2, 0.1, b2, Q0) * exp(TMDEVOL::S(b2, Q0) - TMDEVOL::S(b2, Q[i])));
      }
      p1.SetMarkerColor(colorlist[j]);
      p2.SetMarkerColor(colorlist[j]);
      g1.SetLineColor(colorlist[j]);
      g2.SetLineColor(colorlist[j]);
      c0->cd(1);
      p1.DrawClone("psame");
      g1.DrawClone("lsame");
      c0->cd(2);
      p2.DrawClone("psame");
      g2.DrawClone("lsame");
      j++;
    }
    infile.close();
    
    c0->Print("results/evol_E288_simple.pdf");
  }
  
  
  return 0;
}

