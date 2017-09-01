#include "Lcore.h"
#include "tmdevol.h"
#include "TRandom3.h"

using namespace FIT;

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
    
  return 0;
}
