#include <iostream>

// Has to be placed before the catch2 includes to be registered by them.
std::ostream &operator<<(std::ostream &os, const std::tuple<double, double> &value) {
    os << "azimuth: " << std::get<0>(value) << ", distance: " << std::get<1>(value);
    return os;
}

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#ifdef VINC_C
#include "vinc.h"
#else
#include "vinc.hpp"
#endif

TEST_CASE("Zero distance with zero input") {
    auto [distance, azimuth] = vinc(0, 0, 0, 0);
    REQUIRE(distance == 0);
    REQUIRE(azimuth == 0);
}

TEST_CASE("Zero distance with non zero input") {
    auto [distance, azimuth] = vinc(99, 99, 99, 99);
    REQUIRE(distance == 0);
    REQUIRE(azimuth == 0);
}

TEST_CASE("One degree distance is ~ 111 km") {
    auto [distance, azimuth] = vinc(0, 0, 0, 1);
    REQUIRE_THAT(distance, Catch::Matchers::WithinAbs(111319.491, 0.001));
}

TEST_CASE("Nearly antipodal points") {
    auto [distance, azimuth] = vinc(0, 0.5, 0, 179.5);
    REQUIRE_THAT(distance, Catch::Matchers::WithinAbs(19936288.579, 0.001));
}
