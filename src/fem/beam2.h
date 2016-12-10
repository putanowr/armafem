//! \brief Family of beam2* functions related to 2-dimensional frame element.
//!
//! \copyright GNU General Public License
//! \author Roman Putanowicz

//-----------------------------------------------------------------------------
// This file is part of ArmaFEM project, a C++ library for writing simple
// programs based on Finite Element Method.
// Copyright (C) 2016 Roman Putanowicz <http://foureys.users.greyc.fr/>
//
// ArmaFEM is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ArmaFEM is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// long with ArmaFEM.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------

#pragma once

#include <armadillo>

namespace armafem {

arma::mat beam2e_geom(const arma::vec &n);
arma::mat beam2e_loc(const double E, const double A, const double I, const double L);
void beam2e(const arma::vec &ex, const arma::vec &ey, const arma::vec &ep, arma::mat &Ke);
void beam2e(const arma::vec &ex, const arma::vec &ey, const arma::vec &ep, 
            const arma::vec &eq, arma::mat &Ke, arma::vec &f);

} // namespace armafem
