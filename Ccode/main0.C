#include "Lcore.h"

using namespace FIT;

int main(const int argc, const char * argv[]){

  if (argc < 2) return 0;
  const int opt = atoi(argv[1]);

  SIDIS::FUUT = & SIDIS::Model_FUUT_0;
  PrintLevel = 0;

  //hermes data
  double Q2list[6] = {1.25, 1.51, 1.82, 2.88, 5.23, 9.21};
  SelectionTdelta[1] = 0.1;

  if (opt == 1){
    //Anselmino cut
    SelectionT[0] = 0.5;    SelectionTdelta[0] = 0.5;//x
    SelectionT[2] = 0.3;    SelectionTdelta[2] = 0.3;//z < 0.6
    SelectionT[3] = 0.55;   SelectionTdelta[3] = 0.35;//Pt[0.2,0.9]
  
    FILE * fs = fopen("path/gallery/fitparameters_model0_a.dat","w");
    fprintf(fs, "hermes data proton and deuteron\n");
    fprintf(fs, "z < 0.6, 0.2 < Pt < 0.9\n\n");
    fprintf(fs, "Q2\tNpoints\tChi2\tp0\tp1\n");
    
    double par[2] = {0.5, 0.5};
    for (int i = 0; i < 6; i++){
      Npt = 0;
      SelectionT[1] = Q2list[i];
      
      LoadData("path/Data/SIDIS/hermes.proton.zxpt-3D.vmsub.mults_piplus.list", 19, "proton", "pi+");
      LoadData("path/Data/SIDIS/hermes.proton.zxpt-3D.vmsub.mults_piminus.list", 19, "proton", "pi-");
      LoadData("path/Data/SIDIS/hermes.deuteron.zxpt-3D.vmsub.mults_piplus.list", 19, "deuteron", "pi+");
      LoadData("path/Data/SIDIS/hermes.deuteron.zxpt-3D.vmsub.mults_piminus.list", 19, "deuteron", "pi-");
      
      Minimize(2, par);
      fprintf(fs, "%.2f\t%d\t%.4f\t%.4f\t%.4f\n", SelectionT[1], Npt, Parameters[2], Parameters[0], Parameters[1]);
    }
    fclose(fs);
  }

  if (opt == 2){
    //Anselmino cut 2
    SelectionT[0] = 0.5;    SelectionTdelta[0] = 0.5;//x
    SelectionT[2] = 0.35;    SelectionTdelta[2] = 0.35;//z < 0.7
    SelectionT[3] = 0.55;   SelectionTdelta[3] = 0.35;//Pt[0.2,0.9]
  
    FILE * fs = fopen("path/gallery/fitparameters_model0_b.dat","w");
    fprintf(fs, "hermes data proton and deuteron\n");
    fprintf(fs, "z < 0.7, 0.2 < Pt < 0.9\n\n");
    fprintf(fs, "Q2\tNpoints\tChi2\tp0\tp1\n");
    
    double par[2] = {0.5, 0.5};
    for (int i = 0; i < 6; i++){
      Npt = 0;
      SelectionT[1] = Q2list[i];
      
      LoadData("path/Data/SIDIS/hermes.proton.zxpt-3D.vmsub.mults_piplus.list", 19, "proton", "pi+");
      LoadData("path/Data/SIDIS/hermes.proton.zxpt-3D.vmsub.mults_piminus.list", 19, "proton", "pi-");
      LoadData("path/Data/SIDIS/hermes.deuteron.zxpt-3D.vmsub.mults_piplus.list", 19, "deuteron", "pi+");
      LoadData("path/Data/SIDIS/hermes.deuteron.zxpt-3D.vmsub.mults_piminus.list", 19, "deuteron", "pi-");
      
      Minimize(2, par);
      fprintf(fs, "%.2f\t%d\t%.4f\t%.4f\t%.4f\n", SelectionT[1], Npt, Parameters[2], Parameters[0], Parameters[1]);
    }
    fclose(fs);
  }

  if (opt == 3){
    //Anselmino cut
    SelectionT[0] = 0.5;    SelectionTdelta[0] = 0.5;//x
    SelectionT[2] = 0.3;    SelectionTdelta[2] = 0.3;//z < 0.7
    SelectionT[3] = 0.55;   SelectionTdelta[3] = 5.0;//Pt
  
    FILE * fs = fopen("path/gallery/fitparameters_model0_c.dat","w");
    fprintf(fs, "hermes data proton and deuteron\n");
    fprintf(fs, "z < 0.7, all Pt \n\n");
    fprintf(fs, "Q2\tNpoints\tChi2\tp0\tp1\n");
    
    double par[2] = {0.5, 0.5};
    for (int i = 0; i < 6; i++){
      Npt = 0;
      SelectionT[1] = Q2list[i];
      
      LoadData("path/Data/SIDIS/hermes.proton.zxpt-3D.vmsub.mults_piplus.list", 19, "proton", "pi+");
      LoadData("path/Data/SIDIS/hermes.proton.zxpt-3D.vmsub.mults_piminus.list", 19, "proton", "pi-");
      LoadData("path/Data/SIDIS/hermes.deuteron.zxpt-3D.vmsub.mults_piplus.list", 19, "deuteron", "pi+");
      LoadData("path/Data/SIDIS/hermes.deuteron.zxpt-3D.vmsub.mults_piminus.list", 19, "deuteron", "pi-");
      
      Minimize(2, par);
      fprintf(fs, "%.2f\t%d\t%.4f\t%.4f\t%.4f\n", SelectionT[1], Npt, Parameters[2], Parameters[0], Parameters[1]);
    }
    fclose(fs);
  }

  if (opt == 4){
    //Anselmino cut
    SelectionT[0] = 0.5;    SelectionTdelta[0] = 0.5;//x
    SelectionT[2] = 0.3;    SelectionTdelta[2] = 1.3;//z
    SelectionT[3] = 0.55;   SelectionTdelta[3] = 5.0;//Pt
  
    FILE * fs = fopen("path/gallery/fitparameters_model0_d.dat","w");
    fprintf(fs, "hermes data proton and deuteron\n");
    fprintf(fs, "all z, all Pt \n\n");
    fprintf(fs, "Q2\tNpoints\tChi2\tp0\tp1\n");
    
    double par[2] = {0.5, 0.5};
    for (int i = 0; i < 6; i++){
      Npt = 0;
      SelectionT[1] = Q2list[i];
      
      LoadData("path/Data/SIDIS/hermes.proton.zxpt-3D.vmsub.mults_piplus.list", 19, "proton", "pi+");
      LoadData("path/Data/SIDIS/hermes.proton.zxpt-3D.vmsub.mults_piminus.list", 19, "proton", "pi-");
      LoadData("path/Data/SIDIS/hermes.deuteron.zxpt-3D.vmsub.mults_piplus.list", 19, "deuteron", "pi+");
      LoadData("path/Data/SIDIS/hermes.deuteron.zxpt-3D.vmsub.mults_piminus.list", 19, "deuteron", "pi-");
      
      Minimize(2, par);
      fprintf(fs, "%.2f\t%d\t%.4f\t%.4f\t%.4f\n", SelectionT[1], Npt, Parameters[2], Parameters[0], Parameters[1]);
    }
    fclose(fs);
  }

  

  return 0;
}
