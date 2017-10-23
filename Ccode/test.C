#include "Lcore.h"
#include "tmdevol3.h"

int main(const int argc, const char * argv[]){

  TMDEVOL::xpdf = LHAPDF::mkPDF("CJ15lo", 0);
  cout << TMDEVOL::xpdf->quarkThreshold(4) << endl;

  return 0;
}
