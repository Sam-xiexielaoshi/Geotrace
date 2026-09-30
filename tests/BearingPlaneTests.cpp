#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "GeoTrace/geometry/BearingPlane.h"
#include "GeoTrace/geometry/Vec3.h"

using namespace geotrace::geometry;
using Catch::Matchers::WithinAbs;

TEST_CASE("Bearing plane normal is calculated correctly")
{
    Vec3 observerPosition{
        1.0,
        0.0,
        0.0};

    Vec3 bearingDirection{
        0.0,
        1.0,
        0.0};

    BearingPlane plane = CreateBearingPlane(
        observerPosition,
        bearingDirection);

    REQUIRE_THAT(
        plane.normal.x,
        WithinAbs(0.0, 1e-12));

    REQUIRE_THAT(
        plane.normal.y,
        WithinAbs(0.0, 1e-12));

    REQUIRE_THAT(
        plane.normal.z,
        WithinAbs(1.0, 1e-12));
}

TEST_CASE("Bearing plane normal is perpendicular to position and direction")
{
    Vec3 observerPosition{
        1.0,
        0.0,
        0.0};

    Vec3 bearingDirection{
        0.0,
        1.0,
        0.0};

    BearingPlane plane = CreateBearingPlane(
        observerPosition,
        bearingDirection);

    REQUIRE_THAT(
        Dot(plane.normal, observerPosition),
        WithinAbs(0.0, 1e-12));

    REQUIRE_THAT(
        Dot(plane.normal, bearingDirection),
        WithinAbs(0.0, 1e-12));
}