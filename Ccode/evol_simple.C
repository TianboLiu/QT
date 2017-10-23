#include "LHAPDF/LHAPDF.h"

#include "tmdevol3.h"

#include "TString.h"
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

int main(const int argc, const char * argv[]){

  if (argc < 2){
    cout << "./compare <dataset> <task>" << endl;
    return 0;
  }

  int task = atoi(argv[1]);
  
  if (task == 0){//evol plot
    TMDEVOL::Initialize();
    F_TMD = & F_TMD_evol;

    double bT = 0.8;

    TVirtualPad * d0;
    TCanvas * c0 = new TCanvas("c0", "", 1600, 1800);
    c0->Divide(2,3);

    TH1D * hB = new TH1D("hB", "", 1, 3.0, 20.0);
    hB->SetStats(0);
    hB->SetMinimum(0);
    hB->SetMaximum(10.0);
    hB->GetXaxis()->SetTitle("Q (GeV)");
    hB->GetXaxis()->CenterTitle(true);
    hB->GetXaxis()->SetTitleSize(0.055);
    hB->GetXaxis()->SetTitleOffset(1.15);
    hB->GetXaxis()->SetLabelSize(0.055);
    hB->GetYaxis()->SetTitle("f(x,b,Q)");
    hB->GetYaxis()->CenterTitle(true);
    hB->GetYaxis()->SetTitleSize(0.055);
    hB->GetYaxis()->SetTitleOffset(1.15);
    hB->GetYaxis()->SetLabelSize(0.055);

    TGraph * p0;
    TGraph * g0;

    double LX[100], LY[100];

    for (int i = 0; i < 100; i++)
      LX[i] = 4.0 + 0.15 * i;

    double x = 0.1;

    int flavor = 2;

    char tmp[300];
    double temp;

    x = 0.1;
    d0 = c0->cd(1);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    hB->SetTitle("x=0.1");
    hB->SetMaximum(10.0);
    hB->DrawClone("");
    if (true){
      ifstream infile288("results/QT_DY_E288_simple_selected.txt");
      infile288.getline(tmp, 300);
      while (infile288 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(20);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(4);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(4);
	g0->DrawClone("lsame");
      }
      infile288.close();
      ifstream infile605("results/QT_DY_E605_simple_selected.txt");
      infile605.getline(tmp, 300);
      while (infile605 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(21);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(2);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(2);
	g0->DrawClone("lsame");
      }
      infile605.close();
      ifstream infile772("results/QT_DY_E772_simple_selected.txt");
      infile772.getline(tmp, 300);
      while (infile772 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(22);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(3);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(3);
	g0->DrawClone("lsame");
      }      
      infile772.close();
    }
    
    x = 0.2;
    d0 = c0->cd(2);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    hB->SetTitle("x=0.2");
    hB->SetMaximum(5.0);
    hB->DrawClone("");
    if (true){
      ifstream infile288("results/QT_DY_E288_simple_selected.txt");
      infile288.getline(tmp, 300);
      while (infile288 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(20);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(4);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(4);
	g0->DrawClone("lsame");
      }
      infile288.close();
      ifstream infile605("results/QT_DY_E605_simple_selected.txt");
      infile605.getline(tmp, 300);
      while (infile605 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(21);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(2);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(2);
	g0->DrawClone("lsame");
      }
      infile605.close();
      ifstream infile772("results/QT_DY_E772_simple_selected.txt");
      infile772.getline(tmp, 300);
      while (infile772 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(22);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(3);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(3);
	g0->DrawClone("lsame");
      }
      infile772.close();
    }

    x = 0.3;
    d0 = c0->cd(3);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    hB->SetTitle("x=0.3");
    hB->SetMaximum(2.5);
    hB->DrawClone("");
    if (true){
      ifstream infile288("results/QT_DY_E288_simple_selected.txt");
      infile288.getline(tmp, 300);
      while (infile288 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(20);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(4);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(4);
	g0->DrawClone("lsame");
      }
      infile288.close();
      ifstream infile605("results/QT_DY_E605_simple_selected.txt");
      infile605.getline(tmp, 300);
      while (infile605 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(21);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(2);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(2);
	g0->DrawClone("lsame");
      }
      infile605.close();
      ifstream infile772("results/QT_DY_E772_simple_selected.txt");
      infile772.getline(tmp, 300);
      while (infile772 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(22);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(3);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(3);
	g0->DrawClone("lsame");
      }
      infile772.close();
    }

    x = 0.4;
    d0 = c0->cd(4);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    hB->SetTitle("x=0.4");
    hB->SetMaximum(1.5);
    hB->DrawClone("");
    if (true){
      ifstream infile288("results/QT_DY_E288_simple_selected.txt");
      infile288.getline(tmp, 300);
      while (infile288 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(20);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(4);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(4);
	g0->DrawClone("lsame");
      }
      infile288.close();
      ifstream infile605("results/QT_DY_E605_simple_selected.txt");
      infile605.getline(tmp, 300);
      while (infile605 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(21);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(2);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(2);
	g0->DrawClone("lsame");
      }
      infile605.close();
      ifstream infile772("results/QT_DY_E772_simple_selected.txt");
      infile772.getline(tmp, 300);
      while (infile772 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(22);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(3);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(3);
	g0->DrawClone("lsame");
      }
      infile772.close();
    }

    x = 0.5;
    d0 = c0->cd(5);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    hB->SetTitle("x=0.5");
    hB->SetMaximum(0.8);
    hB->DrawClone("");
    if (true){
      ifstream infile288("results/QT_DY_E288_simple_selected.txt");
      infile288.getline(tmp, 300);
      while (infile288 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(20);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(4);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(4);
	g0->DrawClone("lsame");
      }
      infile288.close();
      ifstream infile605("results/QT_DY_E605_simple_selected.txt");
      infile605.getline(tmp, 300);
      while (infile605 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(21);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(2);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(2);
	g0->DrawClone("lsame");
      }
      infile605.close();
      ifstream infile772("results/QT_DY_E772_simple_selected.txt");
      infile772.getline(tmp, 300);
      while (infile772 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(22);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(3);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(3);
	g0->DrawClone("lsame");
      }
      infile772.close();
    }

    x = 0.6;
    d0 = c0->cd(6);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    hB->SetTitle("x=0.6");
    hB->SetMaximum(0.35);
    hB->DrawClone("");
    if (true){
      ifstream infile288("results/QT_DY_E288_simple_selected.txt");
      infile288.getline(tmp, 300);
      while (infile288 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(20);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(4);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(4);
	g0->DrawClone("lsame");
      }
      infile288.close();
      ifstream infile605("results/QT_DY_E605_simple_selected.txt");
      infile605.getline(tmp, 300);
      while (infile605 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(21);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(2);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(2);
	g0->DrawClone("lsame");
      }
      infile605.close();
      ifstream infile772("results/QT_DY_E772_simple_selected.txt");
      infile772.getline(tmp, 300);
      while (infile772 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(22);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(3);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(3);
	g0->DrawClone("lsame");
      }
      infile772.close();
    }

    c0->Print("results/plot_evol_simple.pdf");
  }

    if (task == 1){//evol plot
    TMDEVOL::Initialize();
    F_TMD = & F_TMD_evol;

    double bT = 0.8;

    TVirtualPad * d0;
    TCanvas * c0 = new TCanvas("c0", "", 1600, 1800);
    c0->Divide(2,3);

    TH1D * hB = new TH1D("hB", "", 1, 3.0, 93.0);
    hB->SetStats(0);
    hB->SetMinimum(0);
    hB->SetMaximum(10.0);
    hB->GetXaxis()->SetTitle("Q (GeV)");
    hB->GetXaxis()->CenterTitle(true);
    hB->GetXaxis()->SetTitleSize(0.055);
    hB->GetXaxis()->SetTitleOffset(1.15);
    hB->GetXaxis()->SetLabelSize(0.055);
    hB->GetYaxis()->SetTitle("f(x,b,Q)");
    hB->GetYaxis()->CenterTitle(true);
    hB->GetYaxis()->SetTitleSize(0.055);
    hB->GetYaxis()->SetTitleOffset(1.15);
    hB->GetYaxis()->SetLabelSize(0.055);

    TGraph * p0;
    TGraph * g0;

    double LX[100], LY[100];

    for (int i = 0; i < 100; i++)
      LX[i] = 3.0 + 0.9 * i;

    double x = 0.1;

    int flavor = 2;

    char tmp[300];
    double temp;

    x = 0.1;
    d0 = c0->cd(1);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    hB->SetTitle("x=0.1");
    hB->SetMaximum(10.0);
    hB->DrawClone("");
    if (true){
      ifstream infile288("results/QT_DY_E288_simple_selected.txt");
      infile288.getline(tmp, 300);
      while (infile288 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(20);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(4);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(4);
	g0->DrawClone("lsame");
      }
      infile288.close();
      ifstream infile605("results/QT_DY_E605_simple_selected.txt");
      infile605.getline(tmp, 300);
      while (infile605 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(21);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(2);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(2);
	g0->DrawClone("lsame");
      }
      infile605.close();
      ifstream infile772("results/QT_DY_E772_simple_selected.txt");
      infile772.getline(tmp, 300);
      while (infile772 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(22);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(3);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(3);
	g0->DrawClone("lsame");
      }      
      infile772.close();

      Q0 = 91.1872;
      Parameters[0] = 2.25;
      Parameters[1] = 0.894;
      Parameters[2] = 2.14;
      p0 = new TGraph(1);
      p0->SetMarkerStyle(24);
      p0->SetMarkerSize(1.5);
      p0->SetMarkerColor(1);
      p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
      p0->DrawClone("psame");
      for (int i = 0; i < 100; i++)
	LY[i] = F_TMD(flavor, x, bT, LX[i]);
      g0 = new TGraph(100, LX, LY);
      g0->SetLineColor(1);
      g0->DrawClone("lsame");		   
    }
    
    x = 0.2;
    d0 = c0->cd(2);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    hB->SetTitle("x=0.2");
    hB->SetMaximum(5.0);
    hB->DrawClone("");
    if (true){
      ifstream infile288("results/QT_DY_E288_simple_selected.txt");
      infile288.getline(tmp, 300);
      while (infile288 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(20);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(4);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(4);
	g0->DrawClone("lsame");
      }
      infile288.close();
      ifstream infile605("results/QT_DY_E605_simple_selected.txt");
      infile605.getline(tmp, 300);
      while (infile605 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(21);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(2);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(2);
	g0->DrawClone("lsame");
      }
      infile605.close();
      ifstream infile772("results/QT_DY_E772_simple_selected.txt");
      infile772.getline(tmp, 300);
      while (infile772 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(22);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(3);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(3);
	g0->DrawClone("lsame");
      }
      infile772.close();

      Q0 = 91.1872;
      Parameters[0] = 2.25;
      Parameters[1] = 0.894;
      Parameters[2] = 2.14;
      p0 = new TGraph(1);
      p0->SetMarkerStyle(24);
      p0->SetMarkerSize(1.5);
      p0->SetMarkerColor(1);
      p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
      p0->DrawClone("psame");
      for (int i = 0; i < 100; i++)
	LY[i] = F_TMD(flavor, x, bT, LX[i]);
      g0 = new TGraph(100, LX, LY);
      g0->SetLineColor(1);
      g0->DrawClone("lsame");
    }

    x = 0.3;
    d0 = c0->cd(3);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    hB->SetTitle("x=0.3");
    hB->SetMaximum(2.5);
    hB->DrawClone("");
    if (true){
      ifstream infile288("results/QT_DY_E288_simple_selected.txt");
      infile288.getline(tmp, 300);
      while (infile288 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(20);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(4);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(4);
	g0->DrawClone("lsame");
      }
      infile288.close();
      ifstream infile605("results/QT_DY_E605_simple_selected.txt");
      infile605.getline(tmp, 300);
      while (infile605 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(21);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(2);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(2);
	g0->DrawClone("lsame");
      }
      infile605.close();
      ifstream infile772("results/QT_DY_E772_simple_selected.txt");
      infile772.getline(tmp, 300);
      while (infile772 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(22);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(3);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(3);
	g0->DrawClone("lsame");
      }
      infile772.close();

      Q0 = 91.1872;
      Parameters[0] = 2.25;
      Parameters[1] = 0.894;
      Parameters[2] = 2.14;
      p0 = new TGraph(1);
      p0->SetMarkerStyle(24);
      p0->SetMarkerSize(1.5);
      p0->SetMarkerColor(1);
      p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
      p0->DrawClone("psame");
      for (int i = 0; i < 100; i++)
	LY[i] = F_TMD(flavor, x, bT, LX[i]);
      g0 = new TGraph(100, LX, LY);
      g0->SetLineColor(1);
      g0->DrawClone("lsame");
    }

    x = 0.4;
    d0 = c0->cd(4);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    hB->SetTitle("x=0.4");
    hB->SetMaximum(1.5);
    hB->DrawClone("");
    if (true){
      ifstream infile288("results/QT_DY_E288_simple_selected.txt");
      infile288.getline(tmp, 300);
      while (infile288 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(20);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(4);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(4);
	g0->DrawClone("lsame");
      }
      infile288.close();
      ifstream infile605("results/QT_DY_E605_simple_selected.txt");
      infile605.getline(tmp, 300);
      while (infile605 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(21);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(2);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(2);
	g0->DrawClone("lsame");
      }
      infile605.close();
      ifstream infile772("results/QT_DY_E772_simple_selected.txt");
      infile772.getline(tmp, 300);
      while (infile772 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(22);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(3);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(3);
	g0->DrawClone("lsame");
      }
      infile772.close();

      Q0 = 91.1872;
      Parameters[0] = 2.25;
      Parameters[1] = 0.894;
      Parameters[2] = 2.14;
      p0 = new TGraph(1);
      p0->SetMarkerStyle(24);
      p0->SetMarkerSize(1.5);
      p0->SetMarkerColor(1);
      p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
      p0->DrawClone("psame");
      for (int i = 0; i < 100; i++)
	LY[i] = F_TMD(flavor, x, bT, LX[i]);
      g0 = new TGraph(100, LX, LY);
      g0->SetLineColor(1);
      g0->DrawClone("lsame");
    }

    x = 0.5;
    d0 = c0->cd(5);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    hB->SetTitle("x=0.5");
    hB->SetMaximum(0.8);
    hB->DrawClone("");
    if (true){
      ifstream infile288("results/QT_DY_E288_simple_selected.txt");
      infile288.getline(tmp, 300);
      while (infile288 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(20);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(4);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(4);
	g0->DrawClone("lsame");
      }
      infile288.close();
      ifstream infile605("results/QT_DY_E605_simple_selected.txt");
      infile605.getline(tmp, 300);
      while (infile605 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(21);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(2);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(2);
	g0->DrawClone("lsame");
      }
      infile605.close();
      ifstream infile772("results/QT_DY_E772_simple_selected.txt");
      infile772.getline(tmp, 300);
      while (infile772 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(22);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(3);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(3);
	g0->DrawClone("lsame");
      }
      infile772.close();

      Q0 = 91.1872;
      Parameters[0] = 2.25;
      Parameters[1] = 0.894;
      Parameters[2] = 2.14;
      p0 = new TGraph(1);
      p0->SetMarkerStyle(24);
      p0->SetMarkerSize(1.5);
      p0->SetMarkerColor(1);
      p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
      p0->DrawClone("psame");
      for (int i = 0; i < 100; i++)
	LY[i] = F_TMD(flavor, x, bT, LX[i]);
      g0 = new TGraph(100, LX, LY);
      g0->SetLineColor(1);
      g0->DrawClone("lsame");
    }

    x = 0.6;
    d0 = c0->cd(6);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    hB->SetTitle("x=0.6");
    hB->SetMaximum(0.35);
    hB->DrawClone("");
    if (true){
      ifstream infile288("results/QT_DY_E288_simple_selected.txt");
      infile288.getline(tmp, 300);
      while (infile288 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(20);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(4);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(4);
	g0->DrawClone("lsame");
      }
      infile288.close();
      ifstream infile605("results/QT_DY_E605_simple_selected.txt");
      infile605.getline(tmp, 300);
      while (infile605 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(21);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(2);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(2);
	g0->DrawClone("lsame");
      }
      infile605.close();
      ifstream infile772("results/QT_DY_E772_simple_selected.txt");
      infile772.getline(tmp, 300);
      while (infile772 >> Q0 >> temp >> temp >> temp >> Parameters[0] >> Parameters[1] >> Parameters[2] >> temp >> temp >> temp){
	p0 = new TGraph(1);
	p0->SetMarkerStyle(22);
	p0->SetMarkerSize(1.5);
	p0->SetMarkerColor(3);
	p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
	p0->DrawClone("psame");
	for (int i = 0; i < 100; i++)
	  LY[i] = F_TMD(flavor, x, bT, LX[i]);
	g0 = new TGraph(100, LX, LY);
	g0->SetLineColor(3);
	g0->DrawClone("lsame");
      }
      infile772.close();

      Q0 = 91.1872;
      Parameters[0] = 2.25;
      Parameters[1] = 0.894;
      Parameters[2] = 2.14;
      p0 = new TGraph(1);
      p0->SetMarkerStyle(24);
      p0->SetMarkerSize(1.5);
      p0->SetMarkerColor(1);
      p0->SetPoint(0, Q0, F_TMD(flavor, x, bT, Q0));
      p0->DrawClone("psame");
      for (int i = 0; i < 100; i++)
	LY[i] = F_TMD(flavor, x, bT, LX[i]);
      g0 = new TGraph(100, LX, LY);
      g0->SetLineColor(1);
      g0->DrawClone("lsame");
    }

    c0->Print("results/plot_evol_simpleZ.pdf");
  }
  

  return 0;
}

