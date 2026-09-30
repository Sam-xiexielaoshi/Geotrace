#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "GeoTrace/geodesy/GeoCoordinate.h"
#include "GeoTrace/geometry/Vec3.h"

using namespace geotrace::geodesy;
using namespace geotrace::geometry;
using Catch::Matchers::WithinAbs;

TEST_CASE("Degrees are converted to radians correctly")
{
    REQUIRE_THAT(
        DegreesToRadians(0.0),
        WithinAbs(0.0, 1e-12));

    REQUIRE_THAT(
        DegreesToRadians(180.0),
        WithinAbs(3.14159265358979323846, 1e-12));
}

TEST_CASE("Equator and prime meridian convert to ECEF correctly")
{
    GeoCoordinate coordinate{
        0.0,
        0.0};

    Vec3 result = LatLonToECEF(coordinate);

    REQUIRE_THAT(result.x, WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(result.y, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(result.z, WithinAbs(0.0, 1e-12));
}

TEST_CASE("Equator and 90 degrees longitude convert to ECEF correctly")
{
    GeoCoordinate coordinate{
        0.0,
        90.0};

    Vec3 result = LatLonToECEF(coordinate);

    REQUIRE_THAT(result.x, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(result.y, WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(result.z, WithinAbs(0.0, 1e-12));
}

TEST_CASE("North pole converts to ECEF correctly")
{
    GeoCoordinate coordinate{
        90.0,
        0.0};

    Vec3 result = LatLonToECEF(coordinate);

    REQUIRE_THAT(result.x, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(result.y, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(result.z, WithinAbs(1.0, 1e-12));
}

TEST_CASE("Local frame at the equator and prime meridian is correct")
{
    GeoCoordinate coordinate{
        0.0,
        0.0};

    LocalFrame frame = ComputeLocalFrame(coordinate);

    REQUIRE_THAT(frame.north.x, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(frame.north.y, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(frame.north.z, WithinAbs(1.0, 1e-12));

    REQUIRE_THAT(frame.east.x, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(frame.east.y, WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(frame.east.z, WithinAbs(0.0, 1e-12));
}

TEST_CASE("Local frame at the equator and 90 degrees longitude is correct")
{
    GeoCoordinate coordinate{
        0.0,
        90.0};

    LocalFrame frame = ComputeLocalFrame(coordinate);

    REQUIRE_THAT(frame.north.x, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(frame.north.y, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(frame.north.z, WithinAbs(1.0, 1e-12));

    REQUIRE_THAT(frame.east.x, WithinAbs(-1.0, 1e-12));
    REQUIRE_THAT(frame.east.y, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(frame.east.z, WithinAbs(0.0, 1e-12));
}

TEST_CASE("Local frame vectors are orthonormal")
{
    GeoCoordinate coordinate{
        20.0,
        40.0};

    LocalFrame frame = ComputeLocalFrame(coordinate);

    REQUIRE_THAT(
        frame.north.Length(),
        WithinAbs(1.0, 1e-12));

    REQUIRE_THAT(
        frame.east.Length(),
        WithinAbs(1.0, 1e-12));

    REQUIRE_THAT(
        Dot(frame.north, frame.east),
        WithinAbs(0.0, 1e-12));
}

TEST_CASE("Bearing of 0 degrees points North")
{
    GeoCoordinate coordinate{
        0.0,
        0.0};

    Vec3 direction = BearingToDirection(
        coordinate,
        0.0);

    REQUIRE_THAT(direction.x, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(direction.y, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(direction.z, WithinAbs(1.0, 1e-12));
}

TEST_CASE("Bearing of 90 degrees points East")
{
    GeoCoordinate coordinate{
        0.0,
        0.0};

    Vec3 direction = BearingToDirection(
        coordinate,
        90.0);

    REQUIRE_THAT(direction.x, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(direction.y, WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(direction.z, WithinAbs(0.0, 1e-12));
}

TEST_CASE("Bearing of 180 degrees points South")
{
    GeoCoordinate coordinate{
        0.0,
        0.0};

    Vec3 direction = BearingToDirection(
        coordinate,
        180.0);

    REQUIRE_THAT(direction.x, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(direction.y, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(direction.z, WithinAbs(-1.0, 1e-12));
}

TEST_CASE("Bearing of 270 degrees points West")
{
    GeoCoordinate coordinate{
        0.0,
        0.0};

    Vec3 direction = BearingToDirection(
        coordinate,
        270.0);

    REQUIRE_THAT(direction.x, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(direction.y, WithinAbs(-1.0, 1e-12));
    REQUIRE_THAT(direction.z, WithinAbs(0.0, 1e-12));
}

TEST_CASE("Bearing direction is a unit vector")
{
    GeoCoordinate coordinate{
        20.0,
        40.0};

    Vec3 direction = BearingToDirection(
        coordinate,
        137.5);

    REQUIRE_THAT(
        direction.Length(),
        WithinAbs(1.0, 1e-12));
}