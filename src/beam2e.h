#pragma once

#include <armadillo>

namespace armafem {

arma::mat beam2e_geom(const arma::vec &n);
arma::mat beam2e_loc(const double E, const double A, const double I, const double L);
void beam2e(const arma::vec &ex, const arma::vec &ey, const arma::vec &ep, arma::mat &Ke);
void beam2e(const arma::vec &ex, const arma::vec &ey, const arma::vec &ep, 
            const arma::vec &eq, arma::mat &Ke, arma::vec &f);

} // namespace armafem
