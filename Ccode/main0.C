#include "Lcore.h"

using namespace FIT;

int main(const int argc, const char * argv[]){

  if (argc < 1) return 0;
  //const int opt = atoi(argv[1]);

  SIDIS::FUUT = & SIDIS::Model_FUUT_0;
  PrintLevel = 0;

  //hermes data
  double Q2list[6] = {1.25, 1.51, 1.82, 2.88, 5.23, 9.21};
  SelectionTdelta[1] = 0.1;

  double Ptmax = 0.2;
  double Ptmin = 0.1;

  //cut
  SelectionT[0] = 0.5;    SelectionTdelta[0] = 0.5;//x
  SelectionT[2] = 0.5;    SelectionTdelta[2] = 0.5;//z
 
  FILE * fs = fopen("path/gallery/ptcut_model0_a.dat","w");
  fprintf(fs, "hermes data proton and deuteron\n");
  fprintf(fs, "z < %.1f\n", SelectionT[2] + SelectionTdelta[2]);
  fprintf(fs, "Q2\t Ptmin\t Ptmax\t Npoints\t Chi2/dof\t p0\t p1\n");
  
  int preNpt = 0;
  for (int i = 0; i < 6; i++){
    SelectionT[1] = Q2list[i];
    preNpt = 0;
    Ptmax = Ptmin;
    while (Ptmax < 2.0){
      Ptmax += 0.05;

      SelectionT[3] = 0.5 * (Ptmax + Ptmin);
      SelectionTdelta[3] = 0.5 * (Ptmax - Ptmin);

      double par[2] = {0.5, 0.5};
   
      Npt = 0;
         
      LoadData("path/Data/SIDIS/hermes.proton.zxpt-3D.vmsub.mults_piplus.list", 19, "proton", "pi+");
      LoadData("path/Data/SIDIS/hermes.proton.zxpt-3D.vmsub.mults_piminus.list", 19, "proton", "pi-");
      LoadData("path/Data/SIDIS/hermes.deuteron.zxpt-3D.vmsub.mults_piplus.list", 19, "deuteron", "pi+");
      LoadData("path/Data/SIDIS/hermes.deuteron.zxpt-3D.vmsub.mults_piminus.list", 19, "deuteron", "pi-");
   
      if (Npt == preNpt) continue;
   
      preNpt = Npt;

      Minimize(2, par);
      fprintf(fs, "%.2f\t %.2f\t %.2f\t %d\t %.2f\t %.4f\t %.4f\n", SelectionT[1], Ptmin, Ptmax, Npt, Parameters[2] / (Npt - 2), Parameters[0], Parameters[1]);
    }
    fprintf(fs, "\n");
  }

  fclose(fs);

  return 0;
}
