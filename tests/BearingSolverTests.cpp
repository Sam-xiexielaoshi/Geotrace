#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "GeoTrace/solver/BearingSolver.h"
#include "GeoTrace/geodesy/GeoCoordinate.h"
#include "GeoTrace/solver/BearingPlaneBuilder.h"
#include "GeoTrace/geometry/Vec3.h"
#include "GeoTrace/solver/BearingObjective.h"
#include "GeoTrace/solver/BearingOptimizer.h"
#include "GeoTrace/solver/BearingIntersectionBuilder.h"
#include "GeoTrace/solver/BearingCandidateCombiner.h"
#include "GeoTrace/solver/BearingObservation.h"
#include "GeoTrace/solver/BearingJacobian.h"
#include "GeoTrace/solver/BearingGaussNewtonOptimizer.h"

#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <random>
#include <string>
#include <iomanip>

using namespace geotrace::solver;
using namespace geotrace::geodesy;
using Catch::Matchers::WithinAbs;
using namespace geotrace::geometry;

namespace
{
    double AddBearingNoise(
        double bearingDegrees,
        double noiseDegrees)
    {
        double noisyBearing =
            bearingDegrees + noiseDegrees;

        // Keep the bearing within [0, 360).
        noisyBearing = std::fmod(
            noisyBearing,
            360.0);

        if (noisyBearing < 0.0)
        {
            noisyBearing += 360.0;
        }

        return noisyBearing;
    }

    double AngularDistanceDegrees(
        const geotrace::geodesy::GeoCoordinate &a,
        const geotrace::geodesy::GeoCoordinate &b)
    {
        const auto positionA =
            geotrace::geodesy::LatLonToECEF(a);

        const auto positionB =
            geotrace::geodesy::LatLonToECEF(b);

        double alignment =
            geotrace::geometry::Dot(positionA, positionB);

        alignment = std::clamp(alignment, -1.0, 1.0);

        constexpr double PI =
            3.14159265358979323846;

        return std::acos(alignment) * 180.0 / PI;
    }
}

TEST_CASE("Bearing plane builder creates one plane per observation")
{
    const GeoCoordinate target{
        20.0,
        40.0};

    const GeoCoordinate observerA{
        0.0,
        0.0};

    const GeoCoordinate observerB{
        0.0,
        90.0};

    const GeoCoordinate observerC{
        30.0,
        45.0};

    const BearingObservation observationA{
        observerA,
        InitialBearing(observerA, target)};

    const BearingObservation observationB{
        observerB,
        InitialBearing(observerB, target)};

    const BearingObservation observationC{
        observerC,
        InitialBearing(observerC, target)};

    const std::vector<BearingObservation> observations{
        observationA,
        observationB,
        observationC};

    const auto planes =
        BuildBearingPlanes(observations);

    REQUIRE(planes.size() == observations.size());

    for (const auto &plane : planes)
    {
        REQUIRE_THAT(
            plane.normal.Length(),
            WithinAbs(1.0, 1e-12));
    }
}

TEST_CASE("Bearing intersection builder creates all unique plane pairs")
{
    const GeoCoordinate target{
        20.0,
        40.0};

    const std::vector<BearingObservation> observations{
        {{0.0, 0.0},
         InitialBearing({0.0, 0.0}, target)},
        {{0.0, 90.0},
         InitialBearing({0.0, 90.0}, target)},
        {{30.0, 45.0},
         InitialBearing({30.0, 45.0}, target)},
        {{-20.0, 10.0},
         InitialBearing({-20.0, 10.0}, target)}};

    const auto planes =
        BuildBearingPlanes(observations);

    const auto intersections =
        BuildBearingIntersections(
            planes,
            observations);

    // 4 observers produce 4 * 3 / 2 = 6 unique pairs.
    REQUIRE(intersections.size() == 6);

    for (const auto &intersection : intersections)
    {
        REQUIRE_THAT(
            intersection.Length(),
            WithinAbs(1.0, 1e-12));
    }
}

TEST_CASE("Bearing candidate combiner averages target directions")
{
    const Vec3 candidateA{
        1.0,
        0.0,
        0.0};

    const Vec3 candidateB{
        1.0,
        0.0,
        0.0};

    const Vec3 candidateC{
        1.0,
        0.0,
        0.0};

    const std::vector<Vec3> candidates{
        candidateA,
        candidateB,
        candidateC};

    const Vec3 result =
        CombineBearingCandidates(candidates);

    REQUIRE_THAT(
        result.x,
        WithinAbs(1.0, 1e-12));

    REQUIRE_THAT(
        result.y,
        WithinAbs(0.0, 1e-12));

    REQUIRE_THAT(
        result.z,
        WithinAbs(0.0, 1e-12));

    REQUIRE_THAT(
        result.Length(),
        WithinAbs(1.0, 1e-12));
}

TEST_CASE("Bearing candidate combiner rejects cancelling candidates")
{
    const std::vector<Vec3> candidates{
        {1.0, 0.0, 0.0},
        {-1.0, 0.0, 0.0}};

    REQUIRE_THROWS(
        CombineBearingCandidates(candidates));
}

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
    const std::vector<BearingObservation> observations{observationA, observationB, observationC};
    GeoCoordinate result = Solve(observations);

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
    const std::vector<BearingObservation> observations{
        observationA,
        observationB,
        observationC};
    GeoCoordinate result = Solve(observations);

    // The solver should reconstruct the new target.
    REQUIRE_THAT(
        result.latitude,
        WithinAbs(35.0, 1e-10));

    REQUIRE_THAT(
        result.longitude,
        WithinAbs(-70.0, 1e-10));
}

TEST_CASE("Four-observer solver reconstructs known target")
{
    const GeoCoordinate target{
        20.0,
        40.0};

    const std::vector<GeoCoordinate> observerLocations{
        {0.0, 0.0},
        {0.0, 90.0},
        {30.0, 45.0},
        {-20.0, 10.0}};

    std::vector<BearingObservation> observations;

    for (const auto &observer : observerLocations)
    {
        observations.push_back({observer,
                                InitialBearing(observer, target)});
    }

    const GeoCoordinate result =
        Solve(observations);

    REQUIRE_THAT(
        result.latitude,
        WithinAbs(target.latitude, 1e-10));

    REQUIRE_THAT(
        result.longitude,
        WithinAbs(target.longitude, 1e-10));
}

TEST_CASE("Five-observer solver reconstructs known target")
{
    const GeoCoordinate target{
        20.0,
        40.0};

    const std::vector<GeoCoordinate> observerLocations{
        {0.0, 0.0},
        {0.0, 90.0},
        {30.0, 45.0},
        {-20.0, 10.0},
        {10.0, -60.0}};

    std::vector<BearingObservation> observations;

    for (const auto &observer : observerLocations)
    {
        observations.push_back({observer,
                                InitialBearing(observer, target)});
    }

    const GeoCoordinate result =
        Solve(observations);

    REQUIRE_THAT(
        result.latitude,
        WithinAbs(target.latitude, 1e-10));

    REQUIRE_THAT(
        result.longitude,
        WithinAbs(target.longitude, 1e-10));
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

    const std::vector<BearingObservation> observations{
        observationA,
        observationB,
        observationC};
    GeoCoordinate result = Solve(observations);

    // Convert the solver result back into 3D.
    Vec3 resultECEF = LatLonToECEF(result);

    // Every point on our spherical Earth model must
    // lie on the unit sphere.
    REQUIRE_THAT(
        resultECEF.Length(),
        WithinAbs(1.0, 1e-12));
}

TEST_CASE("Bearing noise is applied and wrapped correctly")
{
    REQUIRE(
        AddBearingNoise(45.0, 5.0) == 50.0);

    REQUIRE(
        AddBearingNoise(359.0, 3.0) == 2.0);

    REQUIRE(
        AddBearingNoise(1.0, -3.0) == 358.0);
}

TEST_CASE("Three-observer solver responds to different bearing noise levels")
{
    // Known target.
    const GeoCoordinate target{
        20.0,
        40.0};

    // Observer locations.
    const GeoCoordinate observerA{
        0.0,
        0.0};

    const GeoCoordinate observerB{
        0.0,
        90.0};

    const GeoCoordinate observerC{
        30.0,
        45.0};

    // Generate perfect bearings toward the known target.
    const double bearingA =
        InitialBearing(observerA, target);

    const double bearingB =
        InitialBearing(observerB, target);

    const double bearingC =
        InitialBearing(observerC, target);

    // Noise levels we want to study.
    const std::vector<double> noiseLevels{
        0.0,
        0.1,
        0.5,
        1.0,
        2.0,
        5.0};

    for (const double noise : noiseLevels)
    {
        // Apply the same noise pattern to every level.
        const double noisyBearingA =
            AddBearingNoise(bearingA, noise);

        const double noisyBearingB =
            AddBearingNoise(bearingB, -noise);

        const double noisyBearingC =
            AddBearingNoise(bearingC, noise);

        // Build observations.
        const BearingObservation observationA{
            observerA,
            noisyBearingA};

        const BearingObservation observationB{
            observerB,
            noisyBearingB};

        const BearingObservation observationC{
            observerC,
            noisyBearingC};

        const std::vector<BearingObservation> observations{observationA, observationB, observationC};

        // Run the existing solver.
        const GeoCoordinate result =
            Solve(observations);

        // The solver should produce finite coordinates.
        REQUIRE(std::isfinite(result.latitude));
        REQUIRE(std::isfinite(result.longitude));

        // Calculate coordinate error.
        const double latitudeError =
            result.latitude - target.latitude;

        const double longitudeError =
            result.longitude - target.longitude;

        const double positionError =
            std::sqrt(
                latitudeError * latitudeError +
                longitudeError * longitudeError);

        // Print the experiment result.
        // Convert the estimated target back into a 3D position.
        const auto estimatedPosition =
            LatLonToECEF(result);

        // Calculate the bearing residual for each observer.
        const double residualA =
            BearingAngularResidual(
                LatLonToECEF(observerA),
                BearingToDirection(
                    observerA,
                    noisyBearingA),
                estimatedPosition);

        const double residualB =
            BearingAngularResidual(
                LatLonToECEF(observerB),
                BearingToDirection(
                    observerB,
                    noisyBearingB),
                estimatedPosition);

        const double residualC =
            BearingAngularResidual(
                LatLonToECEF(observerC),
                BearingToDirection(
                    observerC,
                    noisyBearingC),
                estimatedPosition);

        // Convert residuals from radians to degrees for readability.
        const double residualADeg =
            residualA * 180.0 / 3.14159265358979323846;

        const double residualBDeg =
            residualB * 180.0 / 3.14159265358979323846;

        const double residualCDeg =
            residualC * 180.0 / 3.14159265358979323846;

        WARN(
            "Noise = " << noise
                       << " deg"
                       << " | Estimated = ("
                       << result.latitude
                       << ", "
                       << result.longitude
                       << ")"
                       << " | Position error = "
                       << positionError
                       << " deg"
                       << " | Residuals = ("
                       << residualADeg
                       << ", "
                       << residualBDeg
                       << ", "
                       << residualCDeg
                       << ") deg");
    }
}

TEST_CASE("Initial bearing and tangent direction agree")
{
    const GeoCoordinate target{
        20.0,
        40.0};

    const GeoCoordinate observer{
        0.0,
        0.0};

    const double bearing =
        InitialBearing(
            observer,
            target);

    const auto observerPosition =
        LatLonToECEF(observer);

    const auto targetPosition =
        LatLonToECEF(target);

    const auto bearingDirection =
        BearingToDirection(
            observer,
            bearing);

    const auto targetDirection =
        TargetTangentDirection(
            observerPosition,
            targetPosition);

    const double alignment =
        Dot(
            bearingDirection,
            targetDirection);

    const double residual =
        BearingAngularResidual(
            observerPosition,
            bearingDirection,
            targetPosition);

    WARN("Initial bearing = " << bearing);
    WARN("Direction alignment = " << alignment);
    WARN("Residual = " << residual * 180.0 / 3.14159265358979323846);

    REQUIRE(
        std::abs(alignment - 1.0) < 1e-12);

    REQUIRE(
        std::abs(residual) < 1e-12);
}

TEST_CASE("Bearing objective is zero for the true target")
{
    const GeoCoordinate target{
        20.0,
        40.0};

    const GeoCoordinate observerA{
        0.0,
        0.0};

    const GeoCoordinate observerB{
        0.0,
        90.0};

    const GeoCoordinate observerC{
        30.0,
        45.0};

    const BearingObservation observationA{
        observerA,
        InitialBearing(observerA, target)};

    const BearingObservation observationB{
        observerB,
        InitialBearing(observerB, target)};

    const BearingObservation observationC{
        observerC,
        InitialBearing(observerC, target)};

    const std::vector<BearingObservation> observations{
        observationA,
        observationB,
        observationC};

    const auto targetPosition =
        LatLonToECEF(target);

    const double error =
        BearingError(
            targetPosition,
            observations);

    INFO("Bearing objective error: " << error);

    REQUIRE(error < 1e-20);
}

TEST_CASE("Bearing objective is larger for a wrong target")
{
    const GeoCoordinate target{
        20.0,
        40.0};

    const GeoCoordinate wrongTarget{
        25.0,
        45.0};

    const GeoCoordinate observerA{
        0.0,
        0.0};

    const GeoCoordinate observerB{
        0.0,
        90.0};

    const GeoCoordinate observerC{
        30.0,
        45.0};

    const BearingObservation observationA{
        observerA,
        InitialBearing(observerA, target)};

    const BearingObservation observationB{
        observerB,
        InitialBearing(observerB, target)};

    const BearingObservation observationC{
        observerC,
        InitialBearing(observerC, target)};

    const std::vector<BearingObservation> observations{
        observationA,
        observationB,
        observationC};

    const auto wrongPosition =
        LatLonToECEF(wrongTarget);

    const double error =
        BearingError(
            wrongPosition,
            observations);

    INFO("Wrong-target bearing objective error: " << error);

    REQUIRE(error > 0.0);
}

TEST_CASE("Bearing objective increases as target moves away")
{
    const geotrace::geodesy::GeoCoordinate trueTarget{
        20.0,
        40.0};

    const geotrace::solver::BearingObservation observationA{
        {0.0, 0.0},
        geotrace::geodesy::InitialBearing(
            {0.0, 0.0},
            trueTarget)};

    const geotrace::solver::BearingObservation observationB{
        {0.0, 90.0},
        geotrace::geodesy::InitialBearing(
            {0.0, 90.0},
            trueTarget)};

    const geotrace::solver::BearingObservation observationC{
        {30.0, 45.0},
        geotrace::geodesy::InitialBearing(
            {30.0, 45.0},
            trueTarget)};

    const std::vector<geotrace::solver::BearingObservation> observations{
        observationA,
        observationB,
        observationC};

    const auto nearbyTarget =
        geotrace::geodesy::LatLonToECEF(
            {21.0, 41.0});

    const auto furtherTarget =
        geotrace::geodesy::LatLonToECEF(
            {30.0, 50.0});

    const auto veryWrongTarget =
        geotrace::geodesy::LatLonToECEF(
            {-20.0, -40.0});

    const double nearbyError =
        geotrace::solver::BearingError(
            nearbyTarget,
            observations);

    const double furtherError =
        geotrace::solver::BearingError(
            furtherTarget,
            observations);

    const double veryWrongError =
        geotrace::solver::BearingError(
            veryWrongTarget,
            observations);

    REQUIRE(nearbyError > 0.0);
    REQUIRE(furtherError > nearbyError);
    REQUIRE(veryWrongError > furtherError);
}

TEST_CASE("Bearing optimizer improves a noisy initial estimate")
{
    const geotrace::geodesy::GeoCoordinate trueTarget{
        20.0,
        40.0};

    const geotrace::solver::BearingObservation observationA{
        {0.0, 0.0},
        geotrace::geodesy::InitialBearing(
            {0.0, 0.0},
            trueTarget) +
            1.0};

    const geotrace::solver::BearingObservation observationB{
        {0.0, 90.0},
        geotrace::geodesy::InitialBearing(
            {0.0, 90.0},
            trueTarget) -
            1.0};

    const geotrace::solver::BearingObservation observationC{
        {30.0, 45.0},
        geotrace::geodesy::InitialBearing(
            {30.0, 45.0},
            trueTarget) +
            1.0};

    const std::vector<geotrace::solver::BearingObservation> observations{
        observationA,
        observationB,
        observationC};

    const auto initialEstimate =
        geotrace::solver::SolveInitial(
            observations);

    const double initialError =
        geotrace::solver::BearingError(
            geotrace::geodesy::LatLonToECEF(
                initialEstimate),
            observations);

    const auto optimizedEstimate =
        geotrace::solver::OptimizeBearingTarget(
            initialEstimate,
            observations);

    const double optimizedError =
        geotrace::solver::BearingError(
            geotrace::geodesy::LatLonToECEF(
                optimizedEstimate),
            observations);

    REQUIRE(optimizedError < initialError);
}

TEST_CASE("Bearing optimizer improves solutions across noise levels")
{
    const geotrace::geodesy::GeoCoordinate trueTarget{
        20.0,
        40.0};

    const std::vector<double> noiseLevels{
        0.0,
        0.1,
        0.5,
        1.0,
        2.0,
        5.0};

    for (const double noise : noiseLevels)
    {
        const geotrace::solver::BearingObservation observationA{
            {0.0, 0.0},
            AddBearingNoise(
                geotrace::geodesy::InitialBearing(
                    {0.0, 0.0},
                    trueTarget),
                noise)};

        const geotrace::solver::BearingObservation observationB{
            {0.0, 90.0},
            AddBearingNoise(
                geotrace::geodesy::InitialBearing(
                    {0.0, 90.0},
                    trueTarget),
                -noise)};

        const geotrace::solver::BearingObservation observationC{
            {30.0, 45.0},
            AddBearingNoise(
                geotrace::geodesy::InitialBearing(
                    {30.0, 45.0},
                    trueTarget),
                noise)};

        const std::vector<geotrace::solver::BearingObservation> observations{
            observationA,
            observationB,
            observationC};

        const auto initialEstimate =
            geotrace::solver::SolveInitial(
                observations);

        const auto optimizedEstimate =
            geotrace::solver::OptimizeBearingTarget(
                initialEstimate,
                observations);

        const double initialError =
            geotrace::solver::BearingError(
                geotrace::geodesy::LatLonToECEF(
                    initialEstimate),
                observations);

        const double optimizedError =
            geotrace::solver::BearingError(
                geotrace::geodesy::LatLonToECEF(
                    optimizedEstimate),
                observations);

        const double initialPositionError =
            AngularDistanceDegrees(
                trueTarget,
                initialEstimate);

        const double optimizedPositionError =
            AngularDistanceDegrees(
                trueTarget,
                optimizedEstimate);

        INFO("Noise = " << noise << " degrees");
        INFO("Initial error = " << initialError);
        INFO("Optimized error = " << optimizedError);

        REQUIRE(optimizedError <= initialError);
        std::cout
            << "Noise = " << noise << " deg\n"
            << "  Initial estimate: ("
            << initialEstimate.latitude << ", "
            << initialEstimate.longitude << ")\n"
            << "  Optimized estimate: ("
            << optimizedEstimate.latitude << ", "
            << optimizedEstimate.longitude << ")\n"
            << "  Initial objective: "
            << initialError << "\n"
            << "  Optimized objective: "
            << optimizedError << "\n"
            << "  Initial position error: "
            << initialPositionError << " deg\n"
            << "  Optimized position error: "
            << optimizedPositionError << " deg\n\n";
    }
}

TEST_CASE("Bearing optimizer is evaluated across randomized noise")
{
    const geotrace::geodesy::GeoCoordinate trueTarget{
        20.0,
        40.0};

    const std::vector<double> noiseLevels{
        0.1,
        0.5,
        1.0,
        2.0,
        5.0};

    constexpr int trialsPerNoiseLevel = 100;

    std::mt19937 generator(42);

    double totalInitialErrorAcrossAllLevels = 0.0;
    double totalOptimizedErrorAcrossAllLevels = 0.0;

    for (const double noiseLevel : noiseLevels)
    {
        std::uniform_real_distribution<double> noiseDistribution(
            -noiseLevel,
            noiseLevel);

        double totalInitialPositionError = 0.0;
        double totalOptimizedPositionError = 0.0;

        int optimizerImprovedCount = 0;

        for (int trial = 0;
             trial < trialsPerNoiseLevel;
             ++trial)
        {
            const geotrace::solver::BearingObservation observationA{
                {0.0, 0.0},
                AddBearingNoise(
                    geotrace::geodesy::InitialBearing(
                        {0.0, 0.0},
                        trueTarget),
                    noiseDistribution(generator))};

            const geotrace::solver::BearingObservation observationB{
                {0.0, 90.0},
                AddBearingNoise(
                    geotrace::geodesy::InitialBearing(
                        {0.0, 90.0},
                        trueTarget),
                    noiseDistribution(generator))};

            const geotrace::solver::BearingObservation observationC{
                {30.0, 45.0},
                AddBearingNoise(
                    geotrace::geodesy::InitialBearing(
                        {30.0, 45.0},
                        trueTarget),
                    noiseDistribution(generator))};

            const std::vector<geotrace::solver::BearingObservation>
                observations{
                    observationA,
                    observationB,
                    observationC};

            const auto initialEstimate =
                geotrace::solver::SolveInitial(
                    observations);

            const auto optimizedEstimate =
                geotrace::solver::OptimizeBearingTarget(
                    initialEstimate,
                    observations);

            const double initialPositionError =
                AngularDistanceDegrees(
                    trueTarget,
                    initialEstimate);

            const double optimizedPositionError =
                AngularDistanceDegrees(
                    trueTarget,
                    optimizedEstimate);

            totalInitialPositionError +=
                initialPositionError;

            totalOptimizedPositionError +=
                optimizedPositionError;

            if (optimizedPositionError <
                initialPositionError)
            {
                ++optimizerImprovedCount;
            }
        }

        const double meanInitialPositionError =
            totalInitialPositionError /
            trialsPerNoiseLevel;

        const double meanOptimizedPositionError =
            totalOptimizedPositionError /
            trialsPerNoiseLevel;

        totalInitialErrorAcrossAllLevels +=
            meanInitialPositionError;

        totalOptimizedErrorAcrossAllLevels +=
            meanOptimizedPositionError;

        std::cout
            << "\nNoise level: "
            << noiseLevel
            << " deg\n"
            << "  Mean initial position error: "
            << meanInitialPositionError
            << " deg\n"
            << "  Mean optimized position error: "
            << meanOptimizedPositionError
            << " deg\n"
            << "  Optimizer improved: "
            << optimizerImprovedCount
            << " / "
            << trialsPerNoiseLevel
            << " trials\n";
    }

    REQUIRE(
        totalOptimizedErrorAcrossAllLevels <
        totalInitialErrorAcrossAllLevels);
}

TEST_CASE("M13.1: Inspect bearing objective landscape")
{
    using namespace geotrace::geodesy;
    using namespace geotrace::solver;

    const GeoCoordinate target{20.0, 40.0};

    const std::vector<GeoCoordinate> observerLocations{
        {0.0, 0.0},
        {0.0, 90.0},
        {30.0, 45.0},
        {-20.0, 10.0},
        {10.0, -60.0}};

    std::vector<BearingObservation> observations;

    for (const auto &observer : observerLocations)
    {
        observations.push_back({observer,
                                InitialBearing(observer, target)});
    }

    INFO("Objective landscape around target (20, 40)");
    INFO("Rows are latitude; columns are longitude.");
    INFO("All objective values are squared angular residuals in radians.");

    for (double latitude = 19.0; latitude <= 21.0001; latitude += 0.5)
    {
        std::string row = "lat=" + std::to_string(latitude) + ": ";

        for (double longitude = 39.0; longitude <= 41.0001; longitude += 0.5)
        {
            const GeoCoordinate candidate{latitude, longitude};

            const double error = BearingError(
                LatLonToECEF(candidate),
                observations);

            row += "(" + std::to_string(longitude) + ", " + std::to_string(error) + ") ";
        }

        WARN(row);
    }

    const double targetError = BearingError(
        LatLonToECEF(target),
        observations);

    REQUIRE(targetError < 1e-12);

    const double offsetError = BearingError(
        LatLonToECEF({20.5, 40.0}),
        observations);

    REQUIRE(offsetError > targetError);
}

TEST_CASE("M13.2: Residual vector contains one value per observer")
{
    using namespace geotrace::geodesy;
    using namespace geotrace::solver;

    const GeoCoordinate target{20.0, 40.0};

    const std::vector<GeoCoordinate> observerLocations{
        {0.0, 0.0},
        {0.0, 90.0},
        {30.0, 45.0},
        {-20.0, 10.0},
        {10.0, -60.0}};

    std::vector<BearingObservation> observations;

    for (const auto &observer : observerLocations)
    {
        observations.push_back({observer,
                                InitialBearing(observer, target)});
    }

    const auto residuals =
        BearingResiduals(LatLonToECEF(target), observations);

    REQUIRE(residuals.size() == observations.size());

    for (const double residual : residuals)
    {
        REQUIRE(std::isfinite(residual));
        REQUIRE(std::abs(residual) < 1e-6);
    }
}

TEST_CASE("M13.2: Objective equals squared residual vector norm")
{
    using namespace geotrace::geodesy;
    using namespace geotrace::solver;

    const GeoCoordinate target{20.0, 40.0};

    const std::vector<GeoCoordinate> observerLocations{
        {0.0, 0.0},
        {0.0, 90.0},
        {30.0, 45.0},
        {-20.0, 10.0},
        {10.0, -60.0}};

    std::vector<BearingObservation> observations;

    for (const auto &observer : observerLocations)
    {
        observations.push_back({observer,
                                InitialBearing(observer, target)});
    }

    // Use a candidate away from the exact solution so the residuals
    // are nonzero and we can test the objective calculation.
    const GeoCoordinate candidate{20.2, 40.3};

    const auto residuals =
        BearingResiduals(LatLonToECEF(candidate), observations);

    double expectedError = 0.0;

    for (const double residual : residuals)
    {
        expectedError += residual * residual;
    }

    const double actualError =
        BearingError(LatLonToECEF(candidate), observations);

    REQUIRE_THAT(actualError, Catch::Matchers::WithinAbs(expectedError, 1e-14));
    REQUIRE(actualError > 0.0);
}

TEST_CASE("M13.3: Numerical residual derivatives are finite and stable")
{
    using namespace geotrace::geodesy;
    using namespace geotrace::solver;

    const GeoCoordinate target{20.0, 40.0};

    const std::vector<GeoCoordinate> observerLocations{
        {0.0, 0.0},
        {0.0, 90.0},
        {30.0, 45.0},
        {-20.0, 10.0},
        {10.0, -60.0}};

    std::vector<BearingObservation> observations;

    for (const auto &observer : observerLocations)
    {
        observations.push_back({observer,
                                InitialBearing(observer, target)});
    }

    // Stay away from the exact solution to avoid testing derivatives
    // at zero angular residual.
    const GeoCoordinate candidate{20.2, 40.3};

    const auto latitudeDerivative =
        [&](std::size_t index, double h)
    {
        const auto north = BearingResiduals(
            LatLonToECEF({candidate.latitude + h,
                          candidate.longitude}),
            observations);

        const auto south = BearingResiduals(
            LatLonToECEF({candidate.latitude - h,
                          candidate.longitude}),
            observations);

        return (north[index] - south[index]) / (2.0 * h);
    };

    const auto longitudeDerivative =
        [&](std::size_t index, double h)
    {
        const auto east = BearingResiduals(
            LatLonToECEF({candidate.latitude,
                          candidate.longitude + h}),
            observations);

        const auto west = BearingResiduals(
            LatLonToECEF({candidate.latitude,
                          candidate.longitude - h}),
            observations);

        return (east[index] - west[index]) / (2.0 * h);
    };

    // Compare two finite-difference step sizes.
    const double coarseStep = 1e-3;
    const double fineStep = 5e-4;

    for (std::size_t i = 0; i < observations.size(); ++i)
    {
        const double dLatCoarse =
            latitudeDerivative(i, coarseStep);

        const double dLatFine =
            latitudeDerivative(i, fineStep);

        const double dLonCoarse =
            longitudeDerivative(i, coarseStep);

        const double dLonFine =
            longitudeDerivative(i, fineStep);

        INFO("Observer index: " << i);

        REQUIRE(std::isfinite(dLatCoarse));
        REQUIRE(std::isfinite(dLatFine));
        REQUIRE(std::isfinite(dLonCoarse));
        REQUIRE(std::isfinite(dLonFine));

        // Derivatives should be reasonably stable when the
        // perturbation step is halved.
        REQUIRE_THAT(
            dLatCoarse,
            Catch::Matchers::WithinAbs(dLatFine, 1e-5));

        REQUIRE_THAT(
            dLonCoarse,
            Catch::Matchers::WithinAbs(dLonFine, 1e-5));
    }
}

TEST_CASE("M13.3.1: Inspect derivative convergence across step sizes")
{
    using namespace geotrace::geodesy;
    using namespace geotrace::solver;

    const GeoCoordinate target{20.0, 40.0};

    const std::vector<GeoCoordinate> observerLocations{
        {0.0, 0.0},
        {0.0, 90.0},
        {30.0, 45.0},
        {-20.0, 10.0},
        {10.0, -60.0}};

    std::vector<BearingObservation> observations;

    for (const auto &observer : observerLocations)
    {
        observations.push_back({observer,
                                InitialBearing(observer, target)});
    }

    const std::vector<GeoCoordinate> candidates{
        {20.2, 40.3},      // Off-target
        {20.001, 40.001}}; // Near-target

    const std::vector<double> stepSizes{
        1e-2,
        1e-3,
        1e-4,
        1e-5};

    // Calculate a central finite-difference derivative
    // for one residual with respect to latitude or longitude.
    const auto derivative =
        [&](const GeoCoordinate &candidate,
            std::size_t residualIndex,
            double h,
            bool latitude)
    {
        GeoCoordinate plus = candidate;
        GeoCoordinate minus = candidate;

        if (latitude)
        {
            plus.latitude += h;
            minus.latitude -= h;
        }
        else
        {
            plus.longitude += h;
            minus.longitude -= h;
        }

        const auto plusResiduals =
            BearingResiduals(LatLonToECEF(plus), observations);

        const auto minusResiduals =
            BearingResiduals(LatLonToECEF(minus), observations);

        return (plusResiduals[residualIndex] - minusResiduals[residualIndex]) / (2.0 * h);
    };

    for (const auto &candidate : candidates)
    {
        WARN("Derivative sweep candidate: latitude = "
             << candidate.latitude
             << ", longitude = "
             << candidate.longitude);

        for (const double h : stepSizes)
        {
            for (std::size_t i = 0; i < observations.size(); ++i)
            {
                const double dLat =
                    derivative(candidate, i, h, true);

                const double dLon =
                    derivative(candidate, i, h, false);

                INFO("Candidate: ("
                     << candidate.latitude << ", "
                     << candidate.longitude << "), h = "
                     << h << ", observer = " << i);

                REQUIRE(std::isfinite(dLat));
                REQUIRE(std::isfinite(dLon));

                WARN("h = " << h
                            << " | observer = " << i
                            << " | dResidual/dLat = " << dLat
                            << " | dResidual/dLon = " << dLon);
            }
        }
    }
}

TEST_CASE("M13.3.2: Signed residual changes sign across the true target")
{
    const GeoCoordinate target{20.0, 40.0};
    const GeoCoordinate observer{0.0, 0.0};

    const std::vector<BearingObservation> observations{
        {observer, InitialBearing(observer, target)}};

    const auto residualAt = [&](const GeoCoordinate &candidate)
    {
        const auto residuals = BearingResiduals(
            LatLonToECEF(candidate),
            observations);

        REQUIRE(residuals.size() == 1);
        return residuals[0];
    };

    const double residualAtTarget = residualAt(target);
    const double residualNorth = residualAt({20.001, 40.0});
    const double residualSouth = residualAt({19.999, 40.0});

    INFO("Residual at target: " << residualAtTarget);
    INFO("Residual north: " << residualNorth);
    INFO("Residual south: " << residualSouth);

    REQUIRE(std::abs(residualAtTarget) < 1e-10);
    REQUIRE(std::isfinite(residualNorth));
    REQUIRE(std::isfinite(residualSouth));
    REQUIRE(std::abs(residualNorth) > 1e-10);
    REQUIRE(std::abs(residualSouth) > 1e-10);
    REQUIRE(residualNorth * residualSouth < 0.0);
}

TEST_CASE("M13.4: Five-point derivatives agree with central differences")
{
    const std::vector<GeoCoordinate> observerLocations{
        {0.0, 0.0},
        {0.0, 90.0},
        {30.0, 45.0},
        {-20.0, 10.0},
        {10.0, -60.0}};

    const GeoCoordinate target{20.0, 40.0};

    std::vector<BearingObservation> observations;
    observations.reserve(observerLocations.size());

    for (const auto &observer : observerLocations)
    {
        observations.push_back({observer,
                                InitialBearing(observer, target)});
    }

    const std::vector<GeoCoordinate> candidates{
        {20.2, 40.3},
        {20.001, 40.001}};

    const auto residualAt = [&](const GeoCoordinate &candidate)
    {
        return BearingResiduals(
            LatLonToECEF(candidate),
            observations);
    };

    const double h = 1e-3;

    for (const auto &candidate : candidates)
    {
        for (int observerIndex = 0;
             observerIndex < static_cast<int>(observations.size());
             ++observerIndex)
        {
            for (const bool differentiateLatitude : {true, false})
            {
                const auto offset = [&](double amount)
                {
                    GeoCoordinate point = candidate;

                    if (differentiateLatitude)
                    {
                        point.latitude += amount;
                    }
                    else
                    {
                        point.longitude += amount;
                    }

                    return point;
                };

                const double fm2 =
                    residualAt(offset(-2.0 * h))[observerIndex];
                const double fm1 =
                    residualAt(offset(-h))[observerIndex];
                const double fp1 =
                    residualAt(offset(h))[observerIndex];
                const double fp2 =
                    residualAt(offset(2.0 * h))[observerIndex];

                const double fivePoint =
                    (fm2 - 8.0 * fm1 + 8.0 * fp1 - fp2) / (12.0 * h);

                const double central =
                    (fp1 - fm1) / (2.0 * h);

                INFO("Candidate: ("
                     << candidate.latitude << ", "
                     << candidate.longitude << ")");
                INFO("Observer index: " << observerIndex);
                INFO("Latitude derivative: "
                     << differentiateLatitude);
                INFO("Five-point derivative: " << fivePoint);
                INFO("Central derivative: " << central);

                REQUIRE(std::isfinite(fivePoint));
                REQUIRE(std::isfinite(central));
                REQUIRE(std::abs(fivePoint - central) < 1e-7);
            }
        }
    }
}

TEST_CASE("M13.5: Jacobian dimensions match the observations")
{
    const GeoCoordinate target{20.0, 40.0};

    const std::vector<GeoCoordinate> observers{
        {0.0, 0.0},
        {0.0, 90.0},
        {30.0, 45.0},
        {-20.0, 10.0},
        {10.0, -60.0}};

    std::vector<BearingObservation> observations;

    for (const auto &observer : observers)
    {
        observations.push_back({observer,
                                InitialBearing(observer, target)});
    }

    const auto jacobian =
        ComputeBearingJacobian(target, observations);

    REQUIRE(jacobian.size() == observations.size());

    for (const auto &row : jacobian)
    {
        REQUIRE(row.size() == 2);
        REQUIRE(std::isfinite(row[0]));
        REQUIRE(std::isfinite(row[1]));
    }
}

TEST_CASE("M13.5: Jacobian rejects non-positive step sizes")
{
    const GeoCoordinate target{20.0, 40.0};

    const std::vector<BearingObservation> observations{
        {{0.0, 0.0}, 60.0}};

    REQUIRE_THROWS_AS(
        ComputeBearingJacobian(target, observations, 0.0),
        std::invalid_argument);

    REQUIRE_THROWS_AS(
        ComputeBearingJacobian(target, observations, -1e-3),
        std::invalid_argument);
}

TEST_CASE("M13.5.1: Jacobian entries agree with five-point derivatives")
{
    const GeoCoordinate target{20.2, 40.3};

    const std::vector<GeoCoordinate> observers{
        {0.0, 0.0},
        {0.0, 90.0},
        {30.0, 45.0},
        {-20.0, 10.0},
        {10.0, -60.0}};

    std::vector<BearingObservation> observations;
    observations.reserve(observers.size());

    for (const auto &observer : observers)
    {
        observations.push_back({observer,
                                InitialBearing(observer, {20.0, 40.0})});
    }

    const double h = 1e-3;

    const auto jacobian =
        ComputeBearingJacobian(target, observations, h);

    const auto residualsAt =
        [&](const GeoCoordinate &point)
    {
        return BearingResiduals(
            LatLonToECEF(point),
            observations);
    };

    for (std::size_t row = 0; row < observations.size(); ++row)
    {
        for (int column = 0; column < 2; ++column)
        {
            const auto offset = [&](double amount)
            {
                auto point = target;

                if (column == 0)
                    point.latitude += amount;
                else
                    point.longitude += amount;

                return point;
            };

            const double fm2 = residualsAt(offset(-2.0 * h))[row];
            const double fm1 = residualsAt(offset(-h))[row];
            const double fp1 = residualsAt(offset(h))[row];
            const double fp2 = residualsAt(offset(2.0 * h))[row];

            const double reference =
                (fm2 - 8.0 * fm1 + 8.0 * fp1 - fp2) / (12.0 * h);

            INFO("Observer row: " << row);
            INFO("Coordinate column: " << column);
            INFO("Jacobian value: " << jacobian[row][column]);
            INFO("Five-point reference: " << reference);

            REQUIRE(std::abs(
                        jacobian[row][column] - reference) < 1e-7);
        }
    }
}

// M13.6.1: Recover an exact target from a perturbed initial estimate.
TEST_CASE("M13.6.1: Gauss-Newton recovers an exact target",
          "[solver][m13][gauss-newton]")
{
    using namespace geotrace;

    const geodesy::GeoCoordinate target{20.0, 40.0};

    const std::vector<geodesy::GeoCoordinate> observerLocations{
        {0.0, 0.0},
        {0.0, 90.0},
        {30.0, 45.0},
        {-20.0, 10.0},
        {10.0, -60.0}};

    std::vector<solver::BearingObservation> observations;

    for (const auto &observer : observerLocations)
    {
        observations.push_back({observer,
                                geodesy::InitialBearing(observer, target)});
    }

    const geodesy::GeoCoordinate initialGuess{22.0, 43.0};

    const auto result =
        solver::OptimizeBearingTargetGaussNewton(
            initialGuess, observations);

    INFO("Initial objective: " << result.initialError);
    INFO("Final objective: " << result.finalError);
    INFO("Estimated latitude: " << result.estimate.latitude);
    INFO("Estimated longitude: " << result.estimate.longitude);

    REQUIRE(result.finalError <= result.initialError);
    REQUIRE(result.finalError < 1e-10);
    REQUIRE(result.estimate.latitude == Catch::Approx(target.latitude)
                                            .margin(1e-3));
    REQUIRE(result.estimate.longitude == Catch::Approx(target.longitude)
                                             .margin(1e-3));
}

// M13.6.1: Reject empty observations and invalid options.
TEST_CASE("M13.6.1: Gauss-Newton validates its inputs",
          "[solver][m13][gauss-newton]")
{
    using namespace geotrace;

    const geodesy::GeoCoordinate initialGuess{20.0, 40.0};
    const std::vector<solver::BearingObservation> noObservations;

    REQUIRE_THROWS_AS(
        solver::OptimizeBearingTargetGaussNewton(
            initialGuess, noObservations),
        std::invalid_argument);

    const std::vector<solver::BearingObservation> observations{
        {{0.0, 0.0}, 45.0}};

    solver::BearingGaussNewtonOptions options;
    options.maxIterations = 0;

    REQUIRE_THROWS_AS(
        solver::OptimizeBearingTargetGaussNewton(
            initialGuess, observations, options),
        std::invalid_argument);

    options = {};
    options.initialDamping = 0.0;

    REQUIRE_THROWS_AS(
        solver::OptimizeBearingTargetGaussNewton(
            initialGuess, observations, options),
        std::invalid_argument);
}

TEST_CASE("M13.6.2: Compare Gauss-Newton with coordinate descent",
          "[solver][m13][gauss-newton]")
{
    using namespace geotrace;

    const geodesy::GeoCoordinate target{20.0, 40.0};

    const std::vector<geodesy::GeoCoordinate> observerLocations{
        {0.0, 0.0},
        {0.0, 90.0},
        {30.0, 45.0},
        {-20.0, 10.0},
        {10.0, -60.0}};

    std::vector<solver::BearingObservation> observations;

    for (const auto &observer : observerLocations)
    {
        observations.push_back({observer,
                                geodesy::InitialBearing(observer, target)});
    }

    // Use the same initial estimate for both algorithms.
    const geodesy::GeoCoordinate initialGuess{22.0, 43.0};

    const auto coordinateDescent =
        solver::OptimizeBearingTarget(
            initialGuess, observations);

    const auto gaussNewton =
        solver::OptimizeBearingTargetGaussNewton(
            initialGuess, observations);

    const double coordinateDescentError =
        solver::BearingError(
            geodesy::LatLonToECEF(coordinateDescent),
            observations);

    const double gaussNewtonError =
        solver::BearingError(
            geodesy::LatLonToECEF(gaussNewton.estimate),
            observations);

    const auto positionErrorDegrees =
        [&](const geodesy::GeoCoordinate &estimate)
    {
        const auto estimatedPosition =
            geodesy::LatLonToECEF(estimate);

        const auto truePosition =
            geodesy::LatLonToECEF(target);

        const double cosine = std::clamp(
            geometry::Dot(
                estimatedPosition.Normalize(),
                truePosition.Normalize()),
            -1.0, 1.0);

        return std::acos(cosine) *
               180.0 / std::acos(-1.0);
    };

    CAPTURE("Target: (20, 40)");
    CAPTURE("Initial guess: (22, 43)");

    std::cout << "\n=== M13.6.2 Optimizer Comparison ===\n"
              << "Target: (20, 40)\n"
              << "Initial guess: (22, 43)\n"
              << "Coordinate descent estimate: ("
              << coordinateDescent.latitude << ", "
              << coordinateDescent.longitude << ")\n"
              << "Gauss-Newton estimate: ("
              << gaussNewton.estimate.latitude << ", "
              << gaussNewton.estimate.longitude << ")\n"
              << "Coordinate descent objective: "
              << coordinateDescentError << '\n'
              << "Gauss-Newton objective: "
              << gaussNewtonError << '\n'
              << "Coordinate descent position error (deg): "
              << positionErrorDegrees(coordinateDescent) << '\n'
              << "Gauss-Newton position error (deg): "
              << positionErrorDegrees(gaussNewton.estimate) << '\n'
              << "Gauss-Newton iterations: "
              << gaussNewton.iterations << '\n'
              << "Gauss-Newton converged: "
              << std::boolalpha << gaussNewton.converged
              << '\n';

    // Both methods must improve or preserve the initial objective.
    const double initialError =
        solver::BearingError(
            geodesy::LatLonToECEF(initialGuess),
            observations);

    CHECK(coordinateDescentError <= initialError);
    CHECK(gaussNewtonError <= initialError);

    // Gauss-Newton should recover this well-conditioned exact solution.
    CHECK(gaussNewtonError < 1e-10);
    CHECK(positionErrorDegrees(gaussNewton.estimate) < 1e-3);

    // Diagnostic comparison, not a requirement that one method
    // always beats the other.
    SUCCEED("Both optimizer results were measured on identical inputs.");
}

TEST_CASE("M13.6.3: Compare optimizers with noisy bearings",
          "[solver][m13][gauss-newton]")
{
    using namespace geotrace;

    const geodesy::GeoCoordinate target{20.0, 40.0};

    const std::vector<geodesy::GeoCoordinate> observers{
        {0.0, 0.0},
        {0.0, 90.0},
        {30.0, 45.0},
        {-20.0, 10.0},
        {10.0, -60.0}};

    const std::vector<double> noiseDegrees{
        0.10, -0.15, 0.20, -0.08, 0.12};

    std::vector<solver::BearingObservation> observations;

    for (std::size_t i = 0; i < observers.size(); ++i)
    {
        observations.push_back({observers[i],
                                geodesy::InitialBearing(observers[i], target) + noiseDegrees[i]});
    }

    const geodesy::GeoCoordinate initialGuess{22.0, 43.0};

    const auto coordinateDescent =
        solver::OptimizeBearingTarget(initialGuess, observations);

    const auto gaussNewton =
        solver::OptimizeBearingTargetGaussNewton(
            initialGuess, observations);

    const auto objective = [&](const geodesy::GeoCoordinate &point)
    {
        return solver::BearingError(
            geodesy::LatLonToECEF(point), observations);
    };

    const auto positionErrorDegrees =
        [&](const geodesy::GeoCoordinate &estimate)
    {
        const auto estimated =
            geodesy::LatLonToECEF(estimate).Normalize();

        const auto actual =
            geodesy::LatLonToECEF(target).Normalize();

        const double cosine = std::clamp(
            geometry::Dot(estimated, actual), -1.0, 1.0);

        return std::acos(cosine) * 180.0 / std::acos(-1.0);
    };

    const double initialError = objective(initialGuess);
    const double coordinateError = objective(coordinateDescent);
    const double gaussNewtonError = objective(gaussNewton.estimate);

    std::cout << "\n=== M13.6.3 Noisy Optimizer Comparison ===\n"
              << "Initial objective: " << initialError << '\n'
              << "Coordinate descent objective: "
              << coordinateError << '\n'
              << "Gauss-Newton objective: "
              << gaussNewtonError << '\n'
              << "Coordinate descent position error (deg): "
              << positionErrorDegrees(coordinateDescent) << '\n'
              << "Gauss-Newton position error (deg): "
              << positionErrorDegrees(gaussNewton.estimate) << '\n'
              << "Gauss-Newton iterations: "
              << gaussNewton.iterations << '\n'
              << "Gauss-Newton converged: "
              << std::boolalpha << gaussNewton.converged << '\n';

    CHECK(coordinateError <= initialError);
    CHECK(gaussNewtonError <= initialError);
    CHECK(std::isfinite(positionErrorDegrees(coordinateDescent)));
    CHECK(std::isfinite(positionErrorDegrees(gaussNewton.estimate)));
}

TEST_CASE("M13.7.1: Gauss-Newton handles multiple initial guesses",
          "[solver][m13][gauss-newton]")
{
    using namespace geotrace;

    const geodesy::GeoCoordinate target{20.0, 40.0};

    const std::vector<geodesy::GeoCoordinate> observers{
        {0.0, 0.0},
        {0.0, 90.0},
        {30.0, 45.0},
        {-20.0, 10.0},
        {10.0, -60.0}};

    std::vector<solver::BearingObservation> observations;

    for (const auto &observer : observers)
    {
        observations.push_back({observer,
                                geodesy::InitialBearing(observer, target)});
    }

    const std::vector<geodesy::GeoCoordinate> initialGuesses{
        {22.0, 43.0},
        {25.0, 50.0},
        {15.0, 30.0},
        {30.0, 20.0}};

    for (const auto &initialGuess : initialGuesses)
    {
        const auto result =
            solver::OptimizeBearingTargetGaussNewton(
                initialGuess, observations);

        INFO("Initial guess: ("
             << initialGuess.latitude << ", "
             << initialGuess.longitude << ")");

        INFO("Final estimate: ("
             << result.estimate.latitude << ", "
             << result.estimate.longitude << ")");

        CHECK(result.finalError <= result.initialError);
        CHECK(result.finalError < 1e-10);
        CHECK(std::abs(
                  result.estimate.latitude - target.latitude) < 1e-3);
        CHECK(std::abs(
                  result.estimate.longitude - target.longitude) < 1e-3);
    }
}

TEST_CASE("M13.7.2: Gauss-Newton handles longitude wrapping",
          "[solver][m13][gauss-newton]")
{
    using namespace geotrace;

    const geodesy::GeoCoordinate target{20.0, 179.5};

    const std::vector<geodesy::GeoCoordinate> observers{
        {0.0, 0.0},
        {0.0, 90.0},
        {30.0, 45.0},
        {-20.0, 10.0},
        {10.0, -60.0}};

    std::vector<solver::BearingObservation> observations;

    for (const auto &observer : observers)
    {
        observations.push_back({observer,
                                geodesy::InitialBearing(observer, target)});
    }

    SECTION("Initial longitude above +180 degrees")
    {
        const auto result =
            solver::OptimizeBearingTargetGaussNewton(
                {22.0, 183.0}, observations);

        CHECK(result.estimate.longitude >= -180.0);
        CHECK(result.estimate.longitude <= 180.0);
        CHECK(result.finalError < 1e-10);
        CHECK(std::abs(result.estimate.latitude - target.latitude) < 1e-3);
        CHECK(std::abs(result.estimate.longitude - target.longitude) < 1e-3);
    }

    SECTION("Initial longitude below -180 degrees")
    {
        const auto result =
            solver::OptimizeBearingTargetGaussNewton(
                {22.0, -183.0}, observations);

        CHECK(result.estimate.longitude >= -180.0);
        CHECK(result.estimate.longitude <= 180.0);
        CHECK(result.finalError < 1e-10);
        CHECK(std::abs(result.estimate.latitude - target.latitude) < 1e-3);
        CHECK(std::abs(result.estimate.longitude - target.longitude) < 1e-3);
    }
}

TEST_CASE("M13.7.3: Gauss-Newton handles nearly aligned observers",
          "[solver][m13][gauss-newton]")
{
    using namespace geotrace;

    const geodesy::GeoCoordinate target{20.0, 40.0};

    // These observers are deliberately close together.
    const std::vector<geodesy::GeoCoordinate> observers{
        {19.0, 39.0},
        {19.1, 39.1},
        {19.2, 39.2},
        {19.3, 39.3},
        {19.4, 39.4}};

    std::vector<solver::BearingObservation> observations;

    for (const auto &observer : observers)
    {
        observations.push_back({observer,
                                geodesy::InitialBearing(observer, target)});
    }

    const auto result =
        solver::OptimizeBearingTargetGaussNewton(
            {22.0, 43.0}, observations);

    INFO("Estimate: ("
         << result.estimate.latitude << ", "
         << result.estimate.longitude << ")");

    INFO("Initial objective: " << result.initialError);
    INFO("Final objective: " << result.finalError);
    INFO("Iterations: " << result.iterations);
    INFO("Converged: " << result.converged);
    std::cout << "\n=== M13.7.3 Geometry Diagnostics ===\n"
              << "Estimate: ("
              << result.estimate.latitude << ", "
              << result.estimate.longitude << ")\n"
              << "Initial objective: " << result.initialError << '\n'
              << "Final objective: " << result.finalError << '\n'
              << "Iterations: " << result.iterations << '\n'
              << "Converged: " << std::boolalpha
              << result.converged << '\n';

    CHECK(std::isfinite(result.estimate.latitude));
    CHECK(std::isfinite(result.estimate.longitude));
    CHECK(std::isfinite(result.finalError));
    CHECK(result.finalError <= result.initialError);
}

TEST_CASE("M13.7.4: Gauss-Newton respects its iteration limit",
          "[solver][m13][gauss-newton]")
{
    using namespace geotrace;

    const geodesy::GeoCoordinate target{20.0, 40.0};

    const std::vector<geodesy::GeoCoordinate> observers{
        {0.0, 0.0},
        {0.0, 90.0},
        {30.0, 45.0},
        {-20.0, 10.0},
        {10.0, -60.0}};

    std::vector<solver::BearingObservation> observations;

    for (const auto &observer : observers)
    {
        observations.push_back({observer,
                                geodesy::InitialBearing(observer, target)});
    }

    solver::BearingGaussNewtonOptions options;
    options.maxIterations = 1;

    const auto result =
        solver::OptimizeBearingTargetGaussNewton(
            {22.0, 43.0}, observations, options);

    CHECK(std::isfinite(result.estimate.latitude));
    CHECK(std::isfinite(result.estimate.longitude));
    CHECK(std::isfinite(result.finalError));
    CHECK(result.iterations <= options.maxIterations);
    CHECK(result.finalError <= result.initialError);
}
