#include "Lcore.h"

using namespace FIT;

int main(const int argc, const char * argv[]){

  SIDIS::FUUT = & SIDIS::Model_FUUT_1;

  Npt = 0;

  SelectionT[0] = 0.15;  SelectionTdelta[0] = 0.01;
  SelectionT[1] = 2.9;   SelectionTdelta[1] = 0.1;
  SelectionT[2] = 0.54;  SelectionTdelta[2] = 0.1;
  SelectionT[3] = 1.0;   SelectionTdelta[3] = 10.0;

  LoadData("path/Data/SIDIS/hermes.proton.zxpt-3D.vmsub.mults_piplus.list", 19, "proton", "pi+");

  double par[4] = {0.25, 0.0, 0.25, 0.0};
  cout << Minimize(4, par) << endl;
  cout << Parameters[0] << "  " << Parameters[1] << "  " << Parameters[2] << "  " << Parameters[3] << endl;

  return 0;
}
