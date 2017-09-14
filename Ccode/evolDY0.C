#include "Lcore.h"
#include "tmdevol2.h"

#include "TStyle.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"


using namespace FITDY;

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

  DY::FUU1DY = & DY::Model_FUU1DY_0;
  PrintLevel = 0;

  //E288 E605
  double Qlist[12] = {4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.0, 11.5, 12.5, 13.5, 15.8};

  if (opt == 2){//evolution

    TMDEVOL::Initialize();
    TMDEVOL::F_input = & TMDEVOL::F_input_model0;

    double x = 0.1;
    double bt = 0.2;
    int color[12] = {1, 4, 2, 6, 1, 4, 2, 6, 1, 4, 2, 6};
    int lstyle[12] = {1, 3, 5, 1, 3, 5, 1, 3, 5, 1, 3, 5};
    double ktfit[12] = {0.749, 0.895, 0.958, 0.940, 0.933, 1.30, 1.06, 1.06, 0.994, 1.07, 3.15, 1.00};
    double dktfit[12] = {3.8e-2, 2.7e-2, 2.6e-2, 2.8e-2, 2.2e-2, 5.9e-2, 5.9e-2, 3.4e-1, 3.4e-2, 2.8e-2, 4.4e-1, 4.3e-2};

    int flavor = 2;

    TGraph * gl = new TGraph(20);
    TGraphErrors * gp = new TGraphErrors(1);
  
    TH1D * hB = new TH1D("hB", "", 1, 2.5, 18.0);
    hB->SetStats(0);
    hB->SetMinimum(0.0);
    hB->SetMaximum(13.0);
    hB->GetXaxis()->SetTitle("Q (GeV)");
    hB->GetXaxis()->CenterTitle(true);
    hB->GetXaxis()->SetTitleSize(0.055);
    hB->GetXaxis()->SetTitleOffset(1.15);
    hB->GetXaxis()->SetLabelSize(0.055);
    hB->GetXaxis()->SetNdivisions(6, 5, 0);
    hB->GetYaxis()->SetTitle("f^{u}(x,b_{T})");
    hB->GetYaxis()->CenterTitle(true);
    hB->GetYaxis()->SetTitleSize(0.055);
    hB->GetYaxis()->SetTitleOffset(1.15);
    hB->GetYaxis()->SetLabelSize(0.055);
    hB->GetYaxis()->SetNdivisions(6, 5, 0);

    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);

    hB->DrawClone("axis");
      
    double Q = 3.0;
    for (int i = 0; i < 12; i++){
      TMDEVOL::kt_model0 = ktfit[i];
      TMDEVOL::Q0_model0 = Qlist[i];
      SetGraph(gl, lstyle[i], color[i], 1., 20, color[i], 1.2);
      SetGraph(gp, 1, color[i], 1., 20, color[i], 1.2);
      gp->SetPoint(0, Qlist[i], TMDEVOL::F_output(flavor, x, bt, Qlist[i]));
      gp->SetPointError(0, 0.0, TMDEVOL::F_output(flavor, x, bt, Qlist[i]) * (0.5 * bt * bt * ktfit[i] * dktfit[i]));
      for (int j = 0; j < 20; j++){
	Q = 3.0 + (18.0 - 3.0) * j;
	gl->SetPoint(j, Q, TMDEVOL::F_output(flavor, x, bt, Q));
      }
      gp->DrawClone("pesame");
      gl->DrawClone("lsame");
    }

    TLatex latex;
    latex.SetTextAlign(12);
    latex.SetTextFont(22);
    latex.SetTextSize(0.05);
    latex.DrawLatex(3.2, 11.5, "x = 0.1, b_{T} = 0.1 GeV^{-1}");

    c0->Print("path/gallery/evolDY0_x0.1_bt0.2_u.pdf");

  }
  
  return 0;
}
