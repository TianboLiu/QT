#include "Lcore.h"

using namespace FIT;

int main(const int argc, const char * argv[]){

  SIDIS::FUUT = & SIDIS::Model_FUUT_0;

  Npt = 0;

  SelectionT[0] = 0.15;  SelectionTdelta[0] = 1.01;
  SelectionT[1] = 10.69;   SelectionTdelta[1] = 9.0;
  SelectionT[2] = 0.3;  SelectionTdelta[2] = 0.3;
  SelectionT[3] = 0.55;   SelectionTdelta[3] = 0.35;

  LoadData("path/Data/SIDIS/hermes.proton.zxpt-3D.vmsub.mults_piplus.list", 19, "proton", "pi+");
  LoadData("path/Data/SIDIS/hermes.proton.zxpt-3D.vmsub.mults_piminus.list", 19, "proton", "pi-");
  LoadData("path/Data/SIDIS/hermes.deuteron.zxpt-3D.vmsub.mults_piplus.list", 19, "deuteron", "pi+");
  LoadData("path/Data/SIDIS/hermes.deuteron.zxpt-3D.vmsub.mults_piminus.list", 19, "deuteron", "pi-");

  double par[4] = {0.25, 0.2, 0.25, 0.0};
  cout << Minimize(2, par) << endl;
  cout << Parameters[0] * Parameters[0] << "  " << Parameters[1] * Parameters[1] << endl;

  //Compare();


  return 0;
}
