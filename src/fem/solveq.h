# pragma once

#include <armadillo>

namespace armafem {
//! Solve static FE-equations considering boundary conditions.
//!
//! \param[in] K       global stiffness matrix, dim(K) = nd x nd
//! \praam[in] f       global load vector, dim(f) = nd x 1
//! \param[in] fixId   indices of constrainded DOFs 
//! \param[in] fixVal  constrainded DOF values 
//! \param[out] d      solution including constrained values dim(d) = nd x 1
//! \param[out] Q      reaction force vector dim(Q) = nd x 1
void solveq(const arma::mat &K, const arma::colvec &f, const arma::uvec& fixId, const arma::colvec &fixVal,
            arma::colvec &d, arma::colvec &Q);
} // namespace armafem
