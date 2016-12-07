#include "fem/bar1e.h"
#include "fem/assem.h"
#include "fem/solveq.h"
#include "io/vtk/vtkExporter.h"
#include "io/gnuplot/gnuplotExporter.h"

#include <iostream>
#include <fstream>

using namespace arma;
using namespace armafem;

int main(int arg, char *argv[]) {

  size_t nElem = 15;                    // number of elements
  auto nNodes = nElem+1;                // number of nodes
  double L=1.0;                         // domain span 
  auto x = linspace<vec>(0, L, nNodes); // nodes' coordinates
  auto nDofs = nNodes;                  // total number of DOFs

  auto K = mat(nDofs, nDofs, fill::zeros); // global stiffness matrix
  auto f = colvec(nDofs, fill::zeros);     // global load vector

  auto edof = umat(nElem,3);
  edof(span::all,0) = regspace<ucolvec>(0,nDofs-2);
  edof(span::all,1) = regspace<ucolvec>(0,nDofs-2);
  edof(span::all,2) = regspace<ucolvec>(1,nDofs-1);

  double dx = x(1) - x(0);

  mat Ke;        // local stiffness matrix
  bar1e(1.0/dx, Ke); // calculate element matrix the same for all elements

  assem(edof, K, Ke);

  auto q = 4.0;
  auto duL = 1.0;
  f(span(1,nDofs-2)).fill(q*dx);
  f(0) = f(nDofs-1) = q*dx/2;
  f(nDofs-1) += duL;

  auto uAtZero = 2.0;
  auto fixId = uvec{0};
  auto fixVal = vec{uAtZero};
  auto u = vec{};
  auto Q = vec{};
  solveq(K,f,fixId,fixVal,u,Q);

  // Output three columns: x, u_FEM, u_exact 
  auto results = mat{nNodes, 3, fill::zeros};
  results(span::all, 0) = x;
  results(span::all, 1) = u;
  results(span::all, 2) = x;
  results(span::all, 2).transform([](double x){ return -2*x*x+5*x+2;});
  
  auto vtk = VtkExporter("elliptic.vtk");
  vtk.exportBar1(x, edof, nullptr, nullptr);

  auto gnuplot = GnuplotExporter("elliptic.gpt");
  auto nodalFields = armafem::FieldsTable{{"value", u}};
  gnuplot.exportBar1(x,edof, &nodalFields);
  std::cout << results << "\n";
  return 0;
}

