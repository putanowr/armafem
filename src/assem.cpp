#include <armadillo>
#include <numeric>

namespace armafem {

using namespace arma;

bool assem(const umat edof, mat &K, const mat &Ke, colvec &f, const colvec &fe) {
  auto const &t = edof.cols(1,edof.n_cols); // skip the first column which is elment numbers
  for (auto i = 0; i<edof.n_rows; ++i) {
    const auto &tr = t.row(i);
    K(tr, tr) += Ke;
    f(tr) += fe;
  }
}

} // namespace armafem

