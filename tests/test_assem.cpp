#include <catch.hpp>
#include <assem.h>

using namespace arma;
using namespace armafem;

TEST_CASE("Assembly test", "[assembly]")
{
  auto Ke = imat{{1,2},{3,4}};
  auto K = imat{4,4, fill::zeros};
  auto edofs = urowvec{{0,2,3}};
  auto Kexpected = imat{4,4,fill::zeros};
  assem(edofs, K, Ke);
  Kexpected(2,2) = 1;
  Kexpected(2,3) = 2;
  Kexpected(3,2) = 3;
  Kexpected(3,3) = 4;
  CHECK(true == approx_equal(Kexpected, K, "absdiff", 0.1));
}
