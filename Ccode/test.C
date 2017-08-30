#include "Lcore.h"
#include "tmdevol.h"

int main(const int argc, const char * argv[]){
/*
  double xf[13];
  DIS::zFF(xf, 0.2, 2.0, "pi+");
  for (int i = 0; i < 13; i++)
    cout << xf[i] << endl;

  cout << endl;

  DIS::zFF(xf, 0.2, 2.0, "pi-");
  for (int i = 0; i < 13; i++)
    cout << xf[i] << endl;

  DIS::zFF(xf, 0.2, 2.0, "pi0");
  for (int i = 0; i < 13; i++)
    cout << xf[i] << endl;
*/

  double var[4] = {0.42, 9.21, 0.42, 0.304};
  double par[2] = {0.5, 0.2};

  SIDIS::FUUT = &SIDIS::Model_FUUT_0;
 
  cout << SIDIS::Multiplicity(var, par, "proton", "pi+") << endl;

  TMDEVOL::Initialize();

  cout << TMDEVOL::bstar(0.5) << endl;

  return 0;
}
