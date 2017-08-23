#include "Lcore.h"

int main(const int argc, const char * argv[]){

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

  return 0;
}
