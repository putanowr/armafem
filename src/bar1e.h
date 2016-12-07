//! PURPOSE
//!  Compute element stiffness matrix
//!  for spring (analog) element.
//!
//! INPUT:  ep = [k];       spring stiffness or analog quantity.
//!
//! OUTPUT: Ke :            stiffness matrix, dim(Ke)= 2 x 2

#include <armadillo>

namespace armafem {

//! \brief Compute element stiffness matrix for spring (analog) element.
//!
//! \param[in]  ep spring stiffness or analog quantity
//! \param[out] Ke stiffness matrix, dim(Ke)= 2 x 2
void bar1e(const double ep, arma::mat &Ke);

} // namespacd armafem
