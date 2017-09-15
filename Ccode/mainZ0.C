#include "Lcore.h"

using namespace FITZ;

int main(const int argc, const char * argv[]){

  if (argc < 2) return 0;
  //const int opt = atoi(argv[1]);

  DY::FUU1Z = & DY::Model_FUU1Z_0;
  PrintLevel = 0;

  //CDF D0
  SelectionT[0] = MZ;
  SelectionTdelta[0] = 0.1;

  double QTmax = 20.0;
  double QTmin = 0.0;

  //cut
  SelectionT[1] = 2.0;    SelectionTdelta[1] = 2.0;//QT
  SelectionT[2] = 0.0;    SelectionTdelta[2] = 20.0;//y
  SelectionT[3] = 1.0e8;  SelectionTdelta[3] = 2.0e8;//s
 
  FILE * fs = fopen("path/gallery/fit_Z_0.dat", "w");
  fprintf(fs, "Q\t QTmin\t QTmax\t Npt\t X2/dof\t p0\t e0\n");
  
  double par[2] = {0.7, 1.4};
   
  Npt = 0;
         
  LoadData("path/Data/DY/DY.CDF_RunI.list", 20, "CDF_RunI");
  //LoadData("path/Data/DY/DY.CDF_RunII.list", 20, "CDF_RunII");
  //LoadData("path/Data/DY/DY.D0_RunI.list", 20, "D0_RunI");
  //LoadData("path/Data/DY/DY.D0_RunII.list", 20, "D0_RunII");
  cout << Npt << endl;
  
  
  Minimize(2, par);
  fprintf(fs, "%.2f\t %.2f\t %.2f\t %d\t %.2f\t %.2E\t %.1E\n", MZ, QTmin, QTmax, Npt, Parameters[2] / (Npt - 1), Parameters[0], ParametersError[0]);

  cout << Parameters[2] / (Npt - 1) << endl;
  cout << Parameters[0] << " " << Parameters[1] << endl;
  
  fclose(fs);

  return 0;
}
