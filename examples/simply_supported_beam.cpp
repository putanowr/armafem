#include "beam2e.h"
#include "assem.h"
#include "solveq.h"

#include <iostream>

using namespace arma;
using namespace armafem;

int main(int arg, char *argv[]) {

  size_t nElem = 15;                     // number of elements
  auto nNodes = nElem+1;                 // number of nodes
  double L=10.0;                         // beam length
  auto x = linspace<vec>(0, L, nNodes);  // nodes' coordinates
  auto nDofs = 3*nNodes;                 // total number of DOFs

  auto K = mat(nDofs, nDofs, fill::zeros); // global stiffness matrix
  auto f = colvec(nDofs, fill::zeros);     // global load vector


  auto ex = vec{x(0), x(1)};  // x-coords for element 0
  auto ey = vec{0.0, 0.0};    // y-coords for element 0

  auto E=200e9;               // Young modulus
  auto A=2e-3;                // cross-section area
  auto I=1.6e-5;	            // innertia moment
  auto ep = vec{E, A, I};     // element properties, the same for all elements

  auto eq = vec{0, -10};      // element load vector, the same for all elements

  mat Ke;    // local stiffness matrix
  colvec fe; // local load vector
  beam2e(ex, ey, ep, eq, Ke, fe); // calculate element matrix and load vector
                                  // the same for all elements

  for (size_t i=0; i<nElem; ++i) {
    auto edof = urowvec{i, 3*i,3*i+1, 3*i+2, 3*i+3, 3*i+4, 3*i+5};
    assem(edof, K, Ke, f, fe);
  }

  auto fixId = uvec{0,1, nDofs-2};
  auto fixVal = vec{0.0, 0.0, 0.0};
  auto d = vec{};
  auto Q = vec{};
  solveq(K,f,fixId,fixVal,d,Q);

  // Outpt two columns x-coors and nodal y-displacement
  auto yofx = mat{nNodes, 2, fill::zeros};
  yofx(span::all, 0) = x;
  yofx(span::all, 1) = d(regspace<uvec>(1,3,nDofs));

  std::cout << yofx << "\n";
  return 0;
}

