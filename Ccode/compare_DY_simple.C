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
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TColor.h"

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
  return result / pow(cos(ATan_b_T), 2);
}

double dsigma_DY(const double Q, const double QT, const double y, const double s){
  double par[4] = {Q, QT, y, s};
  //cout << Q << " " << QT << " " << y << " " << s << endl;
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-6);
  ig.SetFunction(&dsigma_DY_integrand, par);
  double result = ig.Integral(0.0, M_PI_2 - 1e-3);
  //double result = ig.Integral(1e-6, 70.0);
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
    else if (Observable[i] == 1){//E605
      dProtonA = 1.0;
      dProtonB = 29.0 / 64.0;
      double dQ = 0.5;
      if (Variable[i][0] > 11.5) dQ = 1.0;
      if (Variable[i][0] > 13.5) dQ = (18.0 - 13.5) / 2.0;
      theory = dsigma_DY(Variable[i][0], Variable[i][1], Variable[i][2], Variable[i][3])
	* pow(Variable[i][0], 2) * 2.0 * log((Variable[i][0] + dQ) / (Variable[i][0] - dQ)) / M_PI;
      sum += pow(theory - Value[i], 2) / (pow(Errors[i][0], 2) + pow(Errors[i][1], 2));
    }
    else if (Observable[i] == 2){//E772
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
  min->SetMaxIterations(100);
  min->SetPrintLevel(0);
  ROOT::Math::Functor f(&Chi2, NPAR);
  min->SetFunction(f);
  //for (int i = 0; i < NPAR; i++){
  //  min->SetVariable(i, "p", init[i], 1.0e-4);
  //}
  min->SetLimitedVariable(0, "N", init[0], 1e-4, 0.0, 2.0);
  if (option == 0) min->SetFixedVariable(0, "N", 1.0);
  min->SetLowerLimitedVariable(1, "alpha", init[1], 1e-4, 1e-12);
  min->SetLowerLimitedVariable(2, "c", init[2], 1e-4, 1e-12);
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
    cout << "./compare <dataset> <task>" << endl;
    return 0;
  }

  int dataset = atoi(argv[1]);
  int task = atoi(argv[2]);
  option = 1;

  if (dataset == 0){
    LoadData_DY("path/Data/DY/DY.E288_200.list", "E288_200");
    LoadData_DY("path/Data/DY/DY.E288_300.list", "E288_300");
    LoadData_DY("path/Data/DY/DY.E288_400.list", "E288_400");
    //LoadData_DY("path/Data/DY/DY.E605.list", "E605");
    //LoadData_DY("path/Data/DY/DY.E772.list", "E772");
  }
  else if (dataset == 1){
    //LoadData_DY("path/Data/DY/DY.E288_200.list", "E288_200");
    //LoadData_DY("path/Data/DY/DY.E288_300.list", "E288_300");
    //LoadData_DY("path/Data/DY/DY.E288_400.list", "E288_400");
    LoadData_DY("path/Data/DY/DY.E605.list", "E605");
    //LoadData_DY("path/Data/DY/DY.E772.list", "E772");
  }
  else if (dataset == 2){
    //LoadData_DY("path/Data/DY/DY.E288_200.list", "E288_200");
    //LoadData_DY("path/Data/DY/DY.E288_300.list", "E288_300");
    //LoadData_DY("path/Data/DY/DY.E288_400.list", "E288_400");
    //LoadData_DY("path/Data/DY/DY.E605.list", "E605");
    LoadData_DY("path/Data/DY/DY.E772.list", "E772");
  }
  
  if (task == 0){//QT scan
    F_TMD = & F_TMD_simple;
    
    double Q = atof(argv[3]);
    SetFlagQ(Q, 1e-3);
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
      if (option == 0) chi2 = Minimize(NPar, Parameters) / (npt - 2);
      else if (option == 1) chi2 = Minimize(NPar, Parameters) / (npt - 3);
      
      printf("%.1f\t%.1f\t%d\t%.2f\t%.2E\t%.2E\t%.2E\t%.1E\t%.1E\t%.1E\n",
	     Q, QT, npt, chi2, Parameters[0], Parameters[1], Parameters[2], ParametersErrors[0], ParametersErrors[1], ParametersErrors[2]);      
    }
  }

  if (task == 1){//QT-chi2 plot
    TH1D * hB = new TH1D("hB", "", 1, 0.0, 10.0);
    hB->SetStats(0);
    hB->SetMinimum(0);
    hB->SetMaximum(10.0);
    hB->GetXaxis()->SetTitle("Q_{T}^{max} (GeV)");
    hB->GetXaxis()->CenterTitle(true);
    hB->GetXaxis()->SetTitleSize(0.055);
    hB->GetXaxis()->SetTitleOffset(1.15);
    hB->GetXaxis()->SetLabelSize(0.055);
    hB->GetYaxis()->SetTitle("#chi^{2}/dof");
    hB->GetYaxis()->CenterTitle(true);
    hB->GetYaxis()->SetTitleSize(0.055);
    hB->GetYaxis()->SetTitleOffset(1.15);
    hB->GetYaxis()->SetLabelSize(0.055);

    int colorlist[10] = {1, 2, 4, 6, 3, 800, 7, 880, 900, 418};
    int stylelist[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    TGraph * g0[10];
    
    TVirtualPad * d0;
    TCanvas * c0 = new TCanvas("c0", "", 800, 1800);
    c0->Divide(1,3);
    d0 = c0->cd(1);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    hB->SetTitle("E288");
    hB->GetXaxis()->SetLimits(0.0, 6.0);
    hB->SetMaximum(5.0);
    hB->DrawClone("");
    
    int npoints = 0;
    char tmp[300];
    double x[30], y[30], temp;
    double Q, Q0;
    TLegend * leg288 = new TLegend(0.8, 0.2, 0.9, 0.9);
    TString leglist288[10] = {"4.5","5.5","6.5","7.5","8.5","9.5","10.5","11.5","12.5","13.5"};
    for (int i = 0; i < 10; i++){
      Q = 4.5 + i;
      npoints = 0;
      ifstream f288("results/QT_DY_E288_simple.txt");
      f288.getline(tmp, 300);
      while (f288 >> Q0 >> x[npoints] >> temp >> y[npoints] >> temp >> temp >>  temp >> temp >> temp >> temp){
	if (Q == Q0) npoints++;
      }
      f288.close();
      g0[i] = new TGraph(npoints, x, y);
      g0[i]->SetLineColor(colorlist[i]);
      g0[i]->SetLineStyle(stylelist[i]);
      g0[i]->DrawClone("lsame");
      leg288->AddEntry(g0[i], leglist288[i], "l");      
    }
    leg288->DrawClone("same");

    d0 = c0->cd(2);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    hB->SetTitle("E605");
    hB->GetXaxis()->SetLimits(0.0, 6.0);
    hB->SetMaximum(5.0);
    hB->DrawClone("");
    double Qlist605[5] = {7.5, 8.5, 11.0, 12.5, 15.8};
    TLegend * leg605 = new TLegend(0.8, 0.2, 0.9, 0.9);
    TString leglist605[5] = {"7.5", "8.5", "11.0", "12.5", "15.8"};
    for (int i = 0; i < 5; i++){
      Q = Qlist605[i];
      npoints = 0;
      ifstream f605("results/QT_DY_E605_simple.txt");
      f605.getline(tmp, 300);
      while (f605 >> Q0 >> x[npoints] >> temp >> y[npoints] >> temp >> temp >>  temp >> temp >> temp >> temp){
	if (Q == Q0) npoints++;
      }
      f605.close();
      g0[i] = new TGraph(npoints, x, y);
      g0[i]->SetLineColor(colorlist[i]);
      g0[i]->SetLineStyle(stylelist[i]);
      g0[i]->DrawClone("lsame");
      leg605->AddEntry(g0[i], leglist605[i], "l");      
    }
    leg605->DrawClone("same");

    d0 = c0->cd(3);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    hB->SetTitle("E772");
    hB->GetXaxis()->SetLimits(0.0, 6.0);
    hB->SetMaximum(10.0);
    hB->DrawClone("");
    double Qlist772[8] = {5.5, 6.5, 7.5, 8.5, 11.5, 12.5, 13.5, 14.5};
    TLegend * leg772 = new TLegend(0.8, 0.2, 0.9, 0.9);
    TString leglist772[8] = {"5.5", "6.5", "7.5", "8.5", "11.5", "12.5", "13.5", "14.5"};
    for (int i = 0; i < 8; i++){
      Q = Qlist772[i];
      npoints = 0;
      ifstream f772("results/QT_DY_E772_simple.txt");
      f772.getline(tmp, 300);
      while (f772 >> Q0 >> x[npoints] >> temp >> y[npoints] >> temp >> temp >>  temp >> temp >> temp >> temp){
	if (Q == Q0) npoints++;
      }
      f772.close();
      g0[i] = new TGraph(npoints, x, y);
      g0[i]->SetLineColor(colorlist[i]);
      g0[i]->SetLineStyle(stylelist[i]);
      g0[i]->DrawClone("lsame");
      leg772->AddEntry(g0[i], leglist772[i], "l");      
    }
    leg772->DrawClone("same");
    
    c0->Print("results/plot_DY_chi2_vs_QT.pdf");
    
  }

  if (task == 2){//data comparison
    F_TMD = & F_TMD_simple;

    TVirtualPad * d0;
    TCanvas * c0 = new TCanvas("c0", "", 1600, 1800);
    c0->Divide(2, 3);

    TH1D * hB = new TH1D("hB", "", 1, 0.0, 10.0);
    hB->SetStats(0);
    hB->SetMinimum(0);
    hB->SetMaximum(10.0);
    hB->GetXaxis()->SetTitle("Q_{T} (GeV)");
    hB->GetXaxis()->CenterTitle(true);
    hB->GetXaxis()->SetTitleSize(0.055);
    hB->GetXaxis()->SetTitleOffset(1.15);
    hB->GetXaxis()->SetLabelSize(0.055);
    hB->GetYaxis()->SetTitle("Ed#sigma/dq^{3} (GeV^{-4})");
    hB->GetYaxis()->CenterTitle(true);
    hB->GetYaxis()->SetTitleSize(0.055);
    hB->GetYaxis()->SetTitleOffset(1.15);
    hB->GetYaxis()->SetLabelSize(0.055);

    double PX[100], PY[100], EY[100];
    double LX[100], LY[100];
    double Q, QT, temp;
    char tmp[300];

    TGraphErrors * p0;
    TGraph * g0;
  
    int colorlist[10] = {1, 2, 4, 6, 3, 800, 7, 880, 900, 418};
    //int stylelist[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    
    d0 = c0->cd(1);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    d0->SetLogy();
    hB->SetTitle("E288 (200 GeV)");
    hB->GetXaxis()->SetLimits(0.0, 3.5);
    hB->SetMinimum(1e-13);
    hB->SetMaximum(1e-7);
    hB->DrawClone("");

    Npt = 0;
    LoadData_DY("path/Data/DY/DY.E288_200.list", "E288_200");
    double Qlist288[10] = {4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5};
    double QTlist288[10] = {1.9, 2.1, 2.3, 4.8, 4.8, 4.6, 1.9, 3.5, 2.9, 2.5};
    dProtonA = 1.0;
    dProtonB = 78.0 / 195.0;
    for (int i = 0; i < 100; i++)
      LX[i] = i * 0.05;

    int ndata = 0;
    for (int i = 0; i < 7; i++){
      ndata = 0;
      Q = Qlist288[i];
      QT = QTlist288[i];
      SetFlagQ(Q, 1e-3);
      SetFlagQT(QT);
      for (int j = 0; j < Npt; j++){
	if (!(FlagQ[j] && FlagQT[j])) continue;
	PX[ndata] = Variable[j][1];
	PY[ndata] = Value[j];
	EY[ndata] = sqrt(pow(Errors[j][0], 2) + pow(Errors[j][1], 2));
	ndata++;
      }
      p0 = new TGraphErrors(ndata, PX, PY, 0, EY);
      p0->SetMarkerStyle(20);
      p0->SetMarkerSize(0.3);
      p0->SetMarkerColor(colorlist[i]);
      p0->SetLineColor(colorlist[i]);
      p0->DrawClone("pesame");
    }

    for (int i = 0; i < 7; i++){
      ifstream fs("results/QT_DY_E288_simple.txt");
      fs.getline(tmp, 300);
      while (fs >> Q >> QT >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	if (Q == Qlist288[i] && QT == QTlist288[i]) break;
      }
      if (Q == 9.5) continue;
      for (int j = 0; j < 100; j++){
	LY[j] = dsigma_DY(Q, LX[j], 0.40, pow(200.0 + Mp, 2) - pow(200.0, 2))
	  * pow(Q, 2) * 2.0 * log((Q + 0.5) / (Q - 0.5)) / M_PI;
      }
      g0 = new TGraph(100, LX, LY);
      g0->SetLineColor(colorlist[i]);
      g0->DrawClone("lsame");
    }

    //
    d0 = c0->cd(3);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    d0->SetLogy();
    hB->SetTitle("E288 (300 GeV)");
    hB->GetXaxis()->SetLimits(0.0, 3.5);
    hB->SetMinimum(1e-13);
    hB->SetMaximum(1e-7);
    hB->DrawClone("");
    
    Npt = 0;
    LoadData_DY("path/Data/DY/DY.E288_300.list", "E288_300");

    for (int i = 0; i < 8; i++){
      ndata = 0;
      Q = Qlist288[i];
      QT = QTlist288[i];
      SetFlagQ(Q, 1e-3);
      SetFlagQT(QT);
      for (int j = 0; j < Npt; j++){
	if (!(FlagQ[j] && FlagQT[j])) continue;
	PX[ndata] = Variable[j][1];
	PY[ndata] = Value[j];
	EY[ndata] = sqrt(pow(Errors[j][0], 2) + pow(Errors[j][1], 2));
	ndata++;
      }
      p0 = new TGraphErrors(ndata, PX, PY, 0, EY);
      p0->SetMarkerStyle(20);
      p0->SetMarkerSize(0.3);
      p0->SetMarkerColor(colorlist[i]);
      p0->SetLineColor(colorlist[i]);
      p0->DrawClone("pesame");
    }

    for (int i = 0; i < 8; i++){
      ifstream fs("results/QT_DY_E288_simple.txt");
      fs.getline(tmp, 300);
      while (fs >> Q >> QT >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	if (Q == Qlist288[i] && QT == QTlist288[i]) break;
      }
      if (Q == 9.5) continue;
      for (int j = 0; j < 100; j++){
	LY[j] = dsigma_DY(Q, LX[j], 0.21, pow(300.0 + Mp, 2) - pow(300.0, 2))
	  * pow(Q, 2) * 2.0 * log((Q + 0.5) / (Q - 0.5)) / M_PI;
      }
      g0 = new TGraph(100, LX, LY);
      g0->SetLineColor(colorlist[i]);
      g0->DrawClone("lsame");
    }

    //
    d0 = c0->cd(5);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    d0->SetLogy();
    hB->SetTitle("E288 (400 GeV)");
    hB->GetXaxis()->SetLimits(0.0, 3.5);
    hB->SetMinimum(1e-13);
    hB->SetMaximum(1e-7);
    hB->DrawClone("");
    
    Npt = 0;
    LoadData_DY("path/Data/DY/DY.E288_400.list", "E288_400");

    for (int i = 1; i < 10; i++){
      ndata = 0;
      Q = Qlist288[i];
      QT = QTlist288[i];
      SetFlagQ(Q, 1e-3);
      SetFlagQT(QT);
      for (int j = 0; j < Npt; j++){
	if (!(FlagQ[j] && FlagQT[j])) continue;
	PX[ndata] = Variable[j][1];
	PY[ndata] = Value[j];
	EY[ndata] = sqrt(pow(Errors[j][0], 2) + pow(Errors[j][1], 2));
	ndata++;
      }
      p0 = new TGraphErrors(ndata, PX, PY, 0, EY);
      p0->SetMarkerStyle(20);
      p0->SetMarkerSize(0.3);
      p0->SetMarkerColor(colorlist[i]);
      p0->SetLineColor(colorlist[i]);
      p0->DrawClone("pesame");
    }

    for (int i = 1; i < 10; i++){
      ifstream fs("results/QT_DY_E288_simple.txt");
      fs.getline(tmp, 300);
      while (fs >> Q >> QT >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	if (Q == Qlist288[i] && QT == QTlist288[i]) break;
      }
      if (Q == 9.5) continue;
      for (int j = 0; j < 100; j++){
	LY[j] = dsigma_DY(Q, LX[j], 0.03, pow(400.0 + Mp, 2) - pow(400.0, 2))
	  * pow(Q, 2) * 2.0 * log((Q + 0.5) / (Q - 0.5)) / M_PI;
      }
      g0 = new TGraph(100, LX, LY);
      g0->SetLineColor(colorlist[i]);
      g0->DrawClone("lsame");
    }


    //
    d0 = c0->cd(2);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    d0->SetLogy();
    hB->SetTitle("E605 (800 GeV)");
    hB->GetXaxis()->SetLimits(0.0, 3.5);
    hB->SetMinimum(1e-13);
    hB->SetMaximum(1e-7);
    hB->DrawClone("");
    
    Npt = 0;
    LoadData_DY("path/Data/DY/DY.E605.list", "E605");
    double Qlist605[5] = {7.5, 8.5, 11.0, 12.5, 15.8};
    double QTlist605[5] = {1.9, 2.3, 2.3, 2.3, 2.1};

    for (int i = 0; i < 5; i++){
      ndata = 0;
      Q = Qlist605[i];
      QT = QTlist605[i];
      SetFlagQ(Q, 1e-3);
      SetFlagQT(QT);
      for (int j = 0; j < Npt; j++){
	if (!(FlagQ[j] && FlagQT[j])) continue;
	PX[ndata] = Variable[j][1];
	PY[ndata] = Value[j];
	EY[ndata] = sqrt(pow(Errors[j][0], 2) + pow(Errors[j][1], 2));
	ndata++;
      }
      p0 = new TGraphErrors(ndata, PX, PY, 0, EY);
      p0->SetMarkerStyle(20);
      p0->SetMarkerSize(0.3);
      p0->SetMarkerColor(colorlist[i]);
      p0->SetLineColor(colorlist[i]);
      p0->DrawClone("pesame");
    }

    for (int i = 0; i < 5; i++){
      ifstream fs("results/QT_DY_E605_simple.txt");
      fs.getline(tmp, 300);
      while (fs >> Q >> QT >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	if (Q == Qlist605[i] && QT == QTlist605[i]) break;
      }
      for (int j = 0; j < 100; j++){
	double dQ = 0.5;
	double s = pow(800.0 + Mp, 2) - pow(800.0, 2);
	if (Q > 11.5) dQ = 1.0;
	if (Q > 13.5) dQ = (18.0 - 13.5) / 2.0;
	LY[j] = dsigma_DY(Q, LX[j], asinh(sqrt(s) / Q * 0.1 / 2.0), s)
	  * pow(Q, 2) * 2.0 * log((Q + dQ) / (Q - dQ)) / M_PI;
      }
      g0 = new TGraph(100, LX, LY);
      g0->SetLineColor(colorlist[i]);
      g0->DrawClone("lsame");
    }

    //
    d0 = c0->cd(4);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    d0->SetLogy();
    hB->SetTitle("E772 (800 GeV)");
    hB->GetXaxis()->SetLimits(0.0, 3.5);
    hB->SetMinimum(1e-13);
    hB->SetMaximum(1e-7);
    hB->DrawClone("");
    
    Npt = 0;
    LoadData_DY("path/Data/DY/DY.E772.list", "E772");
    double Qlist772[8] = {5.5, 6.5, 7.5, 8.5, 11.5, 12.5, 13.5, 14.5};
    double QTlist772[8] = {2.9, 3.9, 3.4, 3.4, 2.7, 2.7, 1.7, 1.7};

    for (int i = 0; i < 8; i++){
      ndata = 0;
      Q = Qlist772[i];
      QT = QTlist772[i];
      SetFlagQ(Q, 1e-3);
      SetFlagQT(QT);
      for (int j = 0; j < Npt; j++){
	if (!(FlagQ[j] && FlagQT[j])) continue;
	PX[ndata] = Variable[j][1];
	PY[ndata] = Value[j];
	EY[ndata] = sqrt(pow(Errors[j][0], 2) + pow(Errors[j][1], 2));
	ndata++;
      }
      p0 = new TGraphErrors(ndata, PX, PY, 0, EY);
      p0->SetMarkerStyle(20);
      p0->SetMarkerSize(0.3);
      p0->SetMarkerColor(colorlist[i]);
      p0->SetLineColor(colorlist[i]);
      p0->DrawClone("pesame");
    }

    for (int i = 0; i < 8; i++){
      continue;
      ifstream fs("results/QT_DY_E772_simple.txt");
      fs.getline(tmp, 300);
      while (fs >> Q >> QT >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	if (Q == Qlist772[i] && QT == QTlist772[i]) break;
      }
      for (int j = 0; j < 100; j++){
	double dQ = 0.5;
	double s = pow(800.0 + Mp, 2) - pow(800.0, 2);
	LY[j] = dsigma_DY(Q, LX[j], asinh(sqrt(s) / Q * 0.2 / 2.0), s)
	  * pow(Q, 2) * 2.0 * log((Q + dQ) / (Q - dQ)) / M_PI;
      }
      g0 = new TGraph(100, LX, LY);
      g0->SetLineColor(colorlist[i]);
      g0->DrawClone("lsame");
    }

    
    c0->Print("results/plot_DY_compare.pdf");

  }
  

  return 0;
}

