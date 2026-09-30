#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "GeoTrace/solver/BearingSolver.h"
#include "GeoTrace/geodesy/GeoCoordinate.h"
#include "GeoTrace/geometry/Vec3.h"

using namespace geotrace::solver;
using namespace geotrace::geodesy;
using Catch::Matchers::WithinAbs;
using namespace geotrace::geometry;

TEST_CASE("Three-observer solver reconstructs known target")
{
    // Known target.
    GeoCoordinate target{
        20.0,
        40.0};

    // Three observers.
    GeoCoordinate observerA{
        0.0,
        0.0};

    GeoCoordinate observerB{
        0.0,
        90.0};

    GeoCoordinate observerC{
        30.0,
        45.0};

    // Generate the bearings that each observer would see
    // when looking toward the known target.
    double bearingA = InitialBearing(
        observerA,
        target);

    double bearingB = InitialBearing(
        observerB,
        target);

    double bearingC = InitialBearing(
        observerC,
        target);

    // Build the three observations.
    BearingObservation observationA{
        observerA,
        bearingA};

    BearingObservation observationB{
        observerB,
        bearingB};

    BearingObservation observationC{
        observerC,
        bearingC};

    // Run the actual solver.
    GeoCoordinate result = Solve(
        observationA,
        observationB,
        observationC);

    // The reconstructed position should match the known target.
    REQUIRE_THAT(
        result.latitude,
        WithinAbs(20.0, 1e-10));

    REQUIRE_THAT(
        result.longitude,
        WithinAbs(40.0, 1e-10));
}

TEST_CASE("Three-observer solver reconstructs a different target")
{
    // Different known target.
    GeoCoordinate target{
        35.0,
        -70.0};

    // Different observer geometry.
    GeoCoordinate observerA{
        10.0,
        -20.0};

    GeoCoordinate observerB{
        -10.0,
        60.0};

    GeoCoordinate observerC{
        50.0,
        10.0};

    // Generate synthetic bearings toward the known target.
    double bearingA = InitialBearing(
        observerA,
        target);

    double bearingB = InitialBearing(
        observerB,
        target);

    double bearingC = InitialBearing(
        observerC,
        target);

    BearingObservation observationA{
        observerA,
        bearingA};

    BearingObservation observationB{
        observerB,
        bearingB};

    BearingObservation observationC{
        observerC,
        bearingC};

    // Run the solver.
    GeoCoordinate result = Solve(
        observationA,
        observationB,
        observationC);

    // The solver should reconstruct the new target.
    REQUIRE_THAT(
        result.latitude,
        WithinAbs(35.0, 1e-10));

    REQUIRE_THAT(
        result.longitude,
        WithinAbs(-70.0, 1e-10));
}

TEST_CASE("Solver result converts back to a unit ECEF vector")
{
    GeoCoordinate target{
        35.0,
        -70.0};

    GeoCoordinate observerA{
        10.0,
        -20.0};

    GeoCoordinate observerB{
        -10.0,
        60.0};

    GeoCoordinate observerC{
        50.0,
        10.0};

    double bearingA = InitialBearing(
        observerA,
        target);

    double bearingB = InitialBearing(
        observerB,
        target);

    double bearingC = InitialBearing(
        observerC,
        target);

    BearingObservation observationA{
        observerA,
        bearingA};

    BearingObservation observationB{
        observerB,
        bearingB};

    BearingObservation observationC{
        observerC,
        bearingC};

    GeoCoordinate result = Solve(
        observationA,
        observationB,
        observationC);

    // Convert the solver result back into 3D.
    Vec3 resultECEF = LatLonToECEF(result);

    // Every point on our spherical Earth model must
    // lie on the unit sphere.
    REQUIRE_THAT(
        resultECEF.Length(),
        WithinAbs(1.0, 1e-12));
}