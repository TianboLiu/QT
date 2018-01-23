#include <iostream>
#include <fstream>

#include "bessel.h"

using namespace std;

int main(const int argc, const char * argv[]){
  cout << "Testing bessel functions from libraries..." << endl;

  cout << BesselJ_GSL(0, 1e7) << endl;
  return 0;

  FILE * f = fopen("results/j0_test.txt", "w");
  fprintf(f, "x\tIntegral\tGSL\tROOT\tSeries\n");
  
  double x = 0.0;
  while (x < 100.0){
    fprintf(f, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	    x,
	    BesselJ_Integral(0, x),
	    BesselJ_GSL(0, x),
	    BesselJ_ROOT(0, x),
	    BesselJ_Series(0, x));
    x = x + 1.0;
  }

  double p = 2.0;
  while (p < 4.0){
    x = pow(10.0, p);
    fprintf(f, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	    x,
	    BesselJ_Integral(0, x),
	    BesselJ_GSL(0, x),
	    BesselJ_ROOT(0, x),
	    BesselJ_Series(0, x));  
    p = p + 0.1;
  }

  fclose(f);

  return 0;
}
