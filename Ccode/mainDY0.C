#include "Lcore.h"

using namespace FITDY;

int main(const int argc, const char * argv[]){

  if (argc < 2) return 0;
  //const int opt = atoi(argv[1]);

  DY::FUU1DY = & DY::Model_FUU1DY_0;
  PrintLevel = 0;

  //E288, E605
  double Qlist[12] = {4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.0, 11.5, 12.5, 13.5, 15.8};
  SelectionTdelta[0] = 0.1;

  double QTmax = 10.0;
  double QTmin = 0.0;

  //cut
  SelectionT[1] = 5.0;    SelectionTdelta[1] = 5.0;//QT
  SelectionT[2] = 0.0;    SelectionTdelta[2] = 2.0;//y
  SelectionT[3] = 1.0e4;  SelectionTdelta[3] = 2.0e4;//s
 
  FILE * fs = fopen("path/gallery/fit_DY_0.dat", "w");
  fprintf(fs, "Q\t QTmin\t QTmax\t Npt\t X2/dof\t p0\t e0\n");
  
  //int preNpt = 0;
  for (int i = 0; i < 12; i++){
    SelectionT[0] = Qlist[i];//Q

    double par[2] = {0.7, 0.4};
   
    Npt = 0;
         
    LoadData("path/Data/DY/DY.E288_200.list", 15, "E288_200");
    LoadData("path/Data/DY/DY.E288_300.list", 15, "E288_300");
    LoadData("path/Data/DY/DY.E288_400.list", 15, "E288_400");
    //LoadData("path/Data/DY/DY.E605.list", 15, "E605_800");
    cout << Npt << endl;
    
    if (Npt == 0) continue;
   
    Minimize(1, par);
    fprintf(fs, "%.2f\t %.2f\t %.2f\t %d\t %.2f\t %.2E\t %.1E\n", Qlist[i], QTmin, QTmax, Npt, Parameters[1] / (Npt - 1), Parameters[0], ParametersError[0]);
  }

  fclose(fs);

  return 0;
}
