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

    TVirtualPad * d0;
    
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
    hB->SetMaximum(20.0);
    hB->SetMinimum(0.0);

    double tmp, ptmax, chi2;
    TLatex latex;
    latex.SetTextFont(22);
    latex.SetTextSize(0.05);
    latex.SetTextAlign(12);
    TGraph * g0;
    int Clist[6] = {1, 4, 2, 6, 7, 9};

    c0->Divide(2,3);
    d0 = c0->cd(1);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
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
    fs1.close();
    latex.DrawLatex(0.05, 27.0, "z<1.0, P_{T}^{min}=0.1");
    

    d0 = c0->cd(2);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    hB->DrawClone("axis");

    ifstream fs2("path/gallery/ptcut_model0_z1.0_ptmin0.2.dat");
    for (int i = 0; i < 3; i++)
      fs2.getline(ltmp, 300);
    int Q2line2[6] = {16, 14, 11, 9, 10, 10};
       
    for (int i = 0; i < 6; i++){
      g0 = new TGraph(Q2line2[i]);
      for (int j = 0; j < Q2line2[i]; j++){
	fs2 >> tmp >> tmp >> ptmax >> tmp >> chi2 >> tmp >> tmp >> tmp >> tmp;
	g0->SetPoint(j, ptmax, chi2);
      }
      fs2.getline(ltmp, 300);
      SetGraph(g0, 1, Clist[i], 2.0, 20, Clist[i], 0.5);
      g0->DrawClone("lpsame");
      g0->Delete();
    }
    fs2.close();
    latex.DrawLatex(0.05, 27.0, "z<1.0, P_{T}^{min}=0.2");


    d0 = c0->cd(3);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    hB->DrawClone("axis");

    ifstream fs3("path/gallery/ptcut_model0_z0.7_ptmin0.1.dat");
    for (int i = 0; i < 3; i++)
      fs3.getline(ltmp, 300);
    int Q2line3[6] = {15, 14, 11, 10, 11, 12};
       
    for (int i = 0; i < 6; i++){
      g0 = new TGraph(Q2line3[i]);
      for (int j = 0; j < Q2line3[i]; j++){
	fs3 >> tmp >> tmp >> ptmax >> tmp >> chi2 >> tmp >> tmp >> tmp >> tmp;
	g0->SetPoint(j, ptmax, chi2);
      }
      fs3.getline(ltmp, 300);
      SetGraph(g0, 1, Clist[i], 2.0, 20, Clist[i], 0.5);
      g0->DrawClone("lpsame");
      g0->Delete();
    }
    fs3.close();
    latex.DrawLatex(0.05, 27.0, "z<0.7, P_{T}^{min}=0.1");

    d0 = c0->cd(4);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    hB->DrawClone("axis");

    ifstream fs4("path/gallery/ptcut_model0_z0.7_ptmin0.2.dat");
    for (int i = 0; i < 3; i++)
      fs4.getline(ltmp, 300);
    int Q2line4[6] = {14, 13, 10, 9, 10, 10};
       
    for (int i = 0; i < 6; i++){
      g0 = new TGraph(Q2line4[i]);
      for (int j = 0; j < Q2line4[i]; j++){
	fs4 >> tmp >> tmp >> ptmax >> tmp >> chi2 >> tmp >> tmp >> tmp >> tmp;
	g0->SetPoint(j, ptmax, chi2);
      }
      fs4.getline(ltmp, 300);
      SetGraph(g0, 1, Clist[i], 2.0, 20, Clist[i], 0.5);
      g0->DrawClone("lpsame");
      g0->Delete();
    }
    fs4.close();
    latex.DrawLatex(0.05, 27.0, "z<0.7, P_{T}^{min}=0.2");

    d0 = c0->cd(5);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    hB->DrawClone("axis");

    ifstream fs5("path/gallery/ptcut_model0_z0.6_ptmin0.1.dat");
    for (int i = 0; i < 3; i++)
      fs5.getline(ltmp, 300);
    int Q2line5[6] = {13, 13, 11, 10, 11, 12};
       
    for (int i = 0; i < 6; i++){
      g0 = new TGraph(Q2line5[i]);
      for (int j = 0; j < Q2line5[i]; j++){
	fs5 >> tmp >> tmp >> ptmax >> tmp >> chi2 >> tmp >> tmp >> tmp >> tmp;
	g0->SetPoint(j, ptmax, chi2);
      }
      fs5.getline(ltmp, 300);
      SetGraph(g0, 1, Clist[i], 2.0, 20, Clist[i], 0.5);
      g0->DrawClone("lpsame");
      g0->Delete();
    }
    fs5.close();
    latex.DrawLatex(0.05, 27.0, "z<0.6, P_{T}^{min}=0.1");
    

    d0 = c0->cd(6);
    d0->SetLeftMargin(0.15);
    d0->SetBottomMargin(0.15);
    hB->DrawClone("axis");

    ifstream fs6("path/gallery/ptcut_model0_z0.6_ptmin0.2.dat");
    for (int i = 0; i < 3; i++)
      fs6.getline(ltmp, 300);
    int Q2line6[6] = {12, 12, 10, 9, 10, 10};
       
    for (int i = 0; i < 6; i++){
      g0 = new TGraph(Q2line6[i]);
      for (int j = 0; j < Q2line6[i]; j++){
	fs6 >> tmp >> tmp >> ptmax >> tmp >> chi2 >> tmp >> tmp >> tmp >> tmp;
	g0->SetPoint(j, ptmax, chi2);
      }
      fs6.getline(ltmp, 300);
      SetGraph(g0, 1, Clist[i], 2.0, 20, Clist[i], 0.5);
      g0->DrawClone("lpsame");
      g0->Delete();
    }
    fs6.close();
    latex.DrawLatex(0.05, 27.0, "z<0.6, P_{T}^{min}=0.2");

    c0->Print("path/gallery/ptcut_model0.pdf");


  }


  return 0;
}
