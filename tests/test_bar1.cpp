#include <catch.hpp>
#include "bar1e.h" 

using namespace arma;
using namespace armafem;

TEST_CASE("Bar 1D element test", "[stiffness]")
{
  mat Ke;
  auto Kexpected = mat{{2.0, -2.0},{-2.0, 2.0}};
  bar1e(2.0, Ke);
  CHECK(true == approx_equal(Kexpected, Ke, "absdiff", 1e-5));
}
