#include "tmdevol3.h"
#include "load.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLIntegrator.h"

using namespace std;

double Parameters[5] = {0.1, 0.1, 0.1, 0.1, 0.1};

double F_NP_model(const int flavor, const double x, const double b_T){
  double value = exp(-Parameters[0] * pow(b_T, Parameters[1]));
  return Parameters[2] * value;
}

