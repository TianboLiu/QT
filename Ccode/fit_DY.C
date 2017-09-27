#include "tmdevol3.h"
#include "load.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLIntegrator.h"

using namespace std;

double Parameters[5] = {0.1, 0.1, 0.1, 0.1, 0.1};

double F_NP_1(const int flavor, const double x, const double b_T){
  return 1.0;
}

double F_NP_N(const int flavor, const double x, const double b_T){
  return Parameters[0];
}



