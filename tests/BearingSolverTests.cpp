#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "GeoTrace/solver/BearingSolver.h"
#include "GeoTrace/geodesy/GeoCoordinate.h"
#include "GeoTrace/geometry/Vec3.h"
#include "GeoTrace/solver/BearingObjective.h"
#include "GeoTrace/solver/BearingOptimizer.h"

#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <random>

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

        // Run the existing solver.
        const GeoCoordinate result =
            Solve(
                observationA,
                observationB,
                observationC);

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
            observationA,
            observationB,
            observationC);

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
                observationA,
                observationB,
                observationC);

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
                    observationA,
                    observationB,
                    observationC);

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