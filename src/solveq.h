# pragma once

#include <armadillo>

namespace armafem {
bool solveq(const arma::mat &K, const arma::colvec &f, const arma::uvec& pdof, const arma::colvec &dp,
            arma::colvec &d, arma::colvec &Q);
} // namespace armafem
