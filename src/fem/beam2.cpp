#include "beam2.h"

using namespace arma;

namespace armafem {

mat beam2e_geom(const vec &n) {
  auto G = mat{6,6, fill::zeros};
  G(0,0) = G(1,1) = G(3,3) = G(4,4) = n(0);
  G(0,1) = G(3,4) = n(1);
  G(1,0) = G(4,3) = -n(1);
  G(2,2) = G(5,5) = 1.0;
  return G;
}

mat beam2e_loc(const double E, const double A, const double I, const double L) {
  const auto L2 = L*L;
  const auto L3 = L2*L;
  const auto EI = E*I;
  const auto EIL = EI/L;
  const auto EIL2 = EI/L2;
  const auto EIL3 = EI/L3;
  auto EAL = E*A/L;

  auto Kle = mat{ {EAL,        0,       0,  EAL,         0,      0 },
                  {0  ,  12*EIL3,  6*EIL2,    0,  -12*EIL3, 6*EIL2 },
                  {0  ,   6*EIL2,  4*EIL ,    0,   -6*EIL2, 2*EIL  },
                  {EAL,        0,       0,  EAL,         0,      0 },
                  {0  , -12*EIL3, -6*EIL2,    0,   12*EIL3,-6*EIL2 },
                  {0  ,   6*EIL2,  2*EIL ,    0,   -6*EIL2, 4*EIL  }
                };
  return Kle;
}

void beam2e_impl(const vec &ex, const vec &ey, const vec &ep, mat &Ke, const vec *eq=nullptr, vec *f=nullptr) {
  auto b = vec{ex(1)-ex(0), ey(1)-ey(0)};
  auto L = norm(b);
  auto n = b/L;
  auto E=ep(0);  
  auto A=ep(1);  
  auto I=ep(2);
 
  auto Kle = beam2e_loc(E,A,I,L); 
  auto G = beam2e_geom(n);

  Ke=G.t()*Kle*G; 
  if (eq != nullptr && f != nullptr) 
  {
     auto qx = (*eq)(0);
     auto qy = (*eq)(1);
     auto fle = colvec{qx/2, qy/2, qy*L/12, qx/2, qy/2, -qy*L/12};
     fle *= L;
     *f = G.t()*fle;
  }
}

void beam2e(const vec &ex, const vec &ey, const vec &ep, mat &Ke) {
  beam2e_impl(ex, ey, ep, Ke, nullptr, nullptr);
}

void beam2e(const vec &ex, const vec &ey, const vec &ep, const vec &eq, mat &Ke, vec &f) {
  beam2e_impl(ex, ey, ep, Ke, &eq, &f);
}

} // namespace armafem
