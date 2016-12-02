#pragma once

#include <armadillo>

namespace armafem {

template <class Matrix, class Vector>
bool assem(const arma::umat edof, Matrix &K, const Matrix &Ke, 
           Vector &f, const Vector &fe) {
  auto const &t = edof.cols(1,edof.n_cols-1); // skip the first column which is elment numbers
  for (auto i = 0; i<edof.n_rows; ++i) {
    const auto &tr = t.row(i);
    K(tr, tr) += Ke;
    f(tr) += fe;
  }
}

template <class Matrix>
bool assem(const arma::umat edof, Matrix &K, const Matrix &Ke) {
  auto const &t = edof.cols(1,edof.n_cols-1); // skip the first column which is elment numbers
  for (auto i = 0; i<edof.n_rows; ++i) {
    const auto &tr = t.row(i);
    K(tr, tr) += Ke;
  }
}

} // namespace armafem
