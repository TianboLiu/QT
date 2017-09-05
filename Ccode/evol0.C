#include "Lcore.h"
#include "tmdevol2.h"

#include "TStyle.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"


using namespace FIT;

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

  SIDIS::FUUT = & SIDIS::Model_FUUT_0;
  PrintLevel = 0;

  //hermes data
  double Q2list[6] = {1.25, 1.51, 1.82, 2.88, 5.23, 9.21};
  double xlist[6] = {0.0375, 0.0611, 0.0957, 0.152, 0.254, 0.410};
  double Ptmaxlist[6] = {0.3, 0.3, 0.9, 0.9, 0.9, 2.0};//from opt0

  SelectionTdelta[1] = 0.1;

  double zmax = 0.6;

  double Ptmax = 0.2;
  double Ptmin = 0.2;

  //cut
  SelectionT[0] = 0.5;    SelectionTdelta[0] = 0.5;//x
  SelectionT[2] = zmax / 2.0;    SelectionTdelta[2] = zmax / 2.0;//z
 
  if (opt == 0){//find Pt cut for each Q2 bin
    const double chi2max = 2.0;
    
    int preNpt = 0;
    for (int i = 0; i < 6; i++){
      SelectionT[1] = Q2list[i];
      preNpt = 0;
      Ptmax = Ptmin;
      while (Ptmax < 2.0){
	Ptmax += 0.01;
	SelectionT[3] = 0.5 * (Ptmax + Ptmin);
	SelectionTdelta[3] = 0.5 * (Ptmax - Ptmin);
	
	Npt = 0;
	LoadData("path/Data/SIDIS/hermes.proton.zxpt-3D.vmsub.mults_piplus.list", 19, "proton", "pi+");
	LoadData("path/Data/SIDIS/hermes.proton.zxpt-3D.vmsub.mults_piminus.list", 19, "proton", "pi-");
	LoadData("path/Data/SIDIS/hermes.deuteron.zxpt-3D.vmsub.mults_piplus.list", 19, "deuteron", "pi+");
	LoadData("path/Data/SIDIS/hermes.deuteron.zxpt-3D.vmsub.mults_piminus.list", 19, "deuteron", "pi-");
   
	if (Npt == preNpt) continue;   
	preNpt = Npt;
	if (Npt < 10) continue;
	
	double par[2] = {0.7, 0.4};
	Minimize(2, par);
	if (Parameters[2] / (Npt - 2) > chi2max){
	  cout << SelectionT[1] << "  " << Npt << "  " << Ptmax << endl;
	  break;
	}
      } 
    }
  }

  if (opt == 1){//Fit with Ptmax cut
    TRandom3 random(0);
    FILE * fs = fopen("path/gallery/fitresult_z0.6_ptmin0.2_model0.dat", "w");
    fprintf(fs, "Q2\t Ptmin\t Ptmax\t Npt\t X2/dof\t p0\t p1\t e0\t e1\n");
    for (int i = 0; i < 6; i++){
      SelectionT[1] = Q2list[i];
      Ptmax = Ptmaxlist[i];
      SelectionT[3] = (Ptmin + Ptmax) / 2.0;
      SelectionTdelta[3] = (Ptmax - Ptmin) / 2.0;
      
      Npt = 0;
      LoadData("path/Data/SIDIS/hermes.proton.zxpt-3D.vmsub.mults_piplus.list", 19, "proton", "pi+");
      LoadData("path/Data/SIDIS/hermes.proton.zxpt-3D.vmsub.mults_piminus.list", 19, "proton", "pi-");
      LoadData("path/Data/SIDIS/hermes.deuteron.zxpt-3D.vmsub.mults_piplus.list", 19, "deuteron", "pi+");
      LoadData("path/Data/SIDIS/hermes.deuteron.zxpt-3D.vmsub.mults_piminus.list", 19, "deuteron", "pi-");

      double par[2] = {random.Uniform(0.2, 1.0), random.Uniform(0.2, 1.0)};
      Minimize(2, par);

      fprintf(fs, "%.2f\t %.2f\t %.2f\t %d\t %.2E\t %.2E\t %.2E\t %.1E\t %.1E\n",
	      SelectionT[1], Ptmin, Ptmax, Npt, Parameters[2] / (Npt - 2), Parameters[0], Parameters[1], ParametersError[0], ParametersError[1]);
    }
    fclose(fs);
  }


  if (opt == 2){//evolution

    TMDEVOL::Initialize();
    TMDEVOL::F_input = & TMDEVOL::F_input_model0;

    double x = 0.1;
    double bt = 0.1;
    int color[4] = {1, 4, 2, 6};

    TGraph * gl[4];
    TGraph * gp[4];
    for (int i = 0; i < 4; i++){
      gl[i] = new TGraph(20);
      gp[i] = new TGraph(1);
      SetGraph(gl[i], 1, color[i], 1., 20, color[i], 1.2);
      SetGraph(gp[i], 1, color[i], 1., 20, color[i], 1.2);
    }

    int flavor = 2;
    double ktfit[4] = {0.72, 0.788, 0.773, 0.658};
    double Q0input[4] = {sqrt(1.82), sqrt(2.88), sqrt(5.24), sqrt(9.21)};
    double Q = 1.5;
    for (int i = 0; i < 4; i++){
      TMDEVOL::kt_model0 = ktfit[i];
      TMDEVOL::Q0_model0 = Q0input[i];
      gp[i]->SetPoint(0, Q0input[i], TMDEVOL::F_output(flavor, x, bt, Q0input[i]));     
      for (int j = 0; j < 20; j++){
	Q = 1.2 + 0.1 * j;
	gl[i]->SetPoint(j, Q, TMDEVOL::F_output(flavor, x, bt, Q));
      }		
    }

    TH1D * hB = new TH1D("hB", "", 1, 1.0, 3.5);
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
    for (int i = 0; i < 4; i++){
      gp[i]->DrawClone("psame");
      gl[i]->DrawClone("lsame");
    }

    TLatex latex;
    latex.SetTextAlign(12);
    latex.SetTextFont(22);
    latex.SetTextSize(0.05);
    latex.DrawLatex(1.2, 11.5, "x = 0.1, b_{T} = 0.1 GeV^{-1}");

    c0->Print("path/gallery/evol0_x0.1_bt0.1_u.pdf");

  }
  
  return 0;
}
