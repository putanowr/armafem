#include <catch.hpp>


TEST_CASE("Smoke test", "[dummy]")
{
  auto expected = 12;
  auto computed = 12;
  CHECK(expected == computed);
}
