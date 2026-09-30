#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "GeoTrace/geometry/PlaneIntersection.h"
#include "GeoTrace/geometry/BearingPlane.h"
#include "GeoTrace/geometry/Vec3.h"

using namespace geotrace::geometry;
using Catch::Matchers::WithinAbs;

TEST_CASE("Two perpendicular planes intersect correctly")
{
    BearingPlane planeA{
        Vec3{
            1.0,
            0.0,
            0.0}};

    BearingPlane planeB{
        Vec3{
            0.0,
            1.0,
            0.0}};

    Vec3 intersection = IntersectPlanes(
        planeA,
        planeB);

    REQUIRE_THAT(
        intersection.x,
        WithinAbs(0.0, 1e-12));

    REQUIRE_THAT(
        intersection.y,
        WithinAbs(0.0, 1e-12));

    REQUIRE_THAT(
        intersection.z,
        WithinAbs(1.0, 1e-12));
}

TEST_CASE("Plane intersection is perpendicular to both plane normals")
{
    BearingPlane planeA{
        Vec3{
            1.0,
            0.0,
            0.0}};

    BearingPlane planeB{
        Vec3{
            0.0,
            1.0,
            0.0}};

    Vec3 intersection = IntersectPlanes(
        planeA,
        planeB);

    REQUIRE_THAT(
        Dot(intersection, planeA.normal),
        WithinAbs(0.0, 1e-12));

    REQUIRE_THAT(
        Dot(intersection, planeB.normal),
        WithinAbs(0.0, 1e-12));
}