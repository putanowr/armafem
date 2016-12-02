#include <catch.hpp>
#include <iostream>

#include "solveq.h"
#include "assem.h"

using namespace arma;
using namespace armafem;

TEST_CASE("solveq test", "[assembly]")
{
  auto Ke = mat{{4,2},{2,4}};
  auto K = mat{4,4, fill::zeros};
  assem(urowvec{0,0,1}, K, Ke);
  assem(urowvec{1,2,3}, K, Ke);

  auto pdof = uvec{0,3};
  auto dp = colvec{1.0, 1.0};
  auto f = colvec{0, 6,6,0};
  auto d = colvec{};
  auto Q = colvec{};
  armafem::solveq(K, f, pdof, dp, d, Q);
  CHECK(true == approx_equal(d, ones<colvec>(4), "absdiff", 1.0e-5));
}
