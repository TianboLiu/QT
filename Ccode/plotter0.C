#include <iostream>
#include <fstream>
#include <cmath>

#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"

using namespace std;

int SetGraph(TGraph * g, int lstyle, int lcolor, double lwidth, int mstyle, int mcolor, double msize){
  g->SetLineStyle(lstyle);
  g->SetLineColor(lcolor);
  g->SetLineWidth(lwidth);
  g->SetMarkerStyle(mstyle);
  g->SetMarkerColor(mcolor);
  g->SetMarkerSize(msize);
  return 0;
}


int main(const int argc, const char * argv[]){

  if (argc < 2) return 0;
  const int opt = atoi(argv[1]);

  TCanvas * c0 = new TCanvas("c0", "", 850, 1100);
  char ltmp[300];

  if (opt == 0){//ptcut plot

    TH1D * hB = new TH1D("hB", "", 1, 0.0, 1.4);
    hB->SetStats(0);
    hB->GetXaxis()->SetTitle("P_{T}^{max} (GeV)");
    hB->GetXaxis()->CenterTitle(true);
    hB->GetXaxis()->SetTitleSize(0.055);
    hB->GetXaxis()->SetTitleOffset(1.15);
    hB->GetXaxis()->SetLabelSize(0.055);
    hB->GetXaxis()->SetNdivisions(6, 5, 0);
    hB->GetYaxis()->SetTitle("#chi^{2}/dof");
    hB->GetYaxis()->CenterTitle(true);
    hB->GetYaxis()->SetTitleSize(0.055);
    hB->GetYaxis()->SetTitleOffset(1.15);
    hB->GetYaxis()->SetLabelSize(0.055);
    hB->GetYaxis()->SetNdivisions(6, 5, 0);
    hB->SetMaximum(30.0);
    hB->SetMinimum(0.0);

    double tmp, ptmax, chi2;
    TGraph * g0;
    int Clist[6] = {1, 4, 2, 5, 6, 7};

    c0->Divide(2,3);
    c0->cd(1);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    hB->DrawClone("axis");

    ifstream fs1("path/gallery/ptcut_model0_z1.0_ptmin0.1.dat");
    for (int i = 0; i < 3; i++)
      fs1.getline(ltmp, 300);
    int Q2line1[6] = {18, 15, 13, 11, 11, 12};
       
    for (int i = 0; i < 6; i++){
      g0 = new TGraph(Q2line1[i]);
      for (int j = 0; j < Q2line1[i]; j++){
	fs1 >> tmp >> tmp >> ptmax >> tmp >> chi2 >> tmp >> tmp >> tmp >> tmp;
	g0->SetPoint(j, ptmax, chi2);
      }
      fs1.getline(ltmp, 300);
      SetGraph(g0, 1, Clist[i], 2.0, 20, Clist[i], 0.5);
      g0->DrawClone("lpsame");
      g0->Delete();
    }

    


    c0->Print("path/gallery/ptcut_model0.pdf");


  }


  return 0;
}
