#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

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

TEST_CASE("Parallel planes cannot be intersected")
{
    const geotrace::geometry::BearingPlane planeA{
        {1.0, 0.0, 0.0}};

    const geotrace::geometry::BearingPlane planeB{
        {1.0, 0.0, 0.0}};

    REQUIRE_THROWS(
        geotrace::geometry::IntersectPlanes(
            planeA,
            planeB));
}

TEST_CASE("Intersection direction resolves the correct antipodal candidate")
{
    using geotrace::geodesy::GeoCoordinate;
    using geotrace::geometry::Vec3;
    using geotrace::solver::BearingObservation;

    const BearingObservation observationA{
        {0.0, 0.0},
        90.0};

    const BearingObservation observationB{
        {0.0, 180.0},
        270.0};

    const Vec3 intersection{
        0.0, 1.0, 0.0};

    const Vec3 resolved =
        geotrace::geometry::ResolveIntersectionDirection(
            intersection,
            observationA,
            observationB);

    REQUIRE(std::abs(resolved.x - 0.0) < 1e-12);
    REQUIRE(std::abs(resolved.y - 1.0) < 1e-12);
    REQUIRE(std::abs(resolved.z - 0.0) < 1e-12);
}