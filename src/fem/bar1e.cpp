#include "bar1e.h"

using namespace arma;

namespace armafem {

void bar1e(const double ep, mat &Ke) {
  Ke = mat{{ep, -ep}, {-ep, ep}}; 
}

} // namespace armafem
