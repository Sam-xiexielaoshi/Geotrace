
#pragma once

#include "GeoTrace/geodesy/GeoCoordinate.h"
#include "GeoTrace/solver/BearingObservation.h"

#include <cstddef>
#include <vector>

namespace geotrace::solver
{
    struct BearingGaussNewtonOptions
    {
        std::size_t maxIterations = 100;
        double initialDamping = 1e-3;
        double gradientTolerance = 1e-10;
        double stepToleranceDegrees = 1e-8;
    };

    struct BearingGaussNewtonResult
    {
        geodesy::GeoCoordinate estimate;
        double initialError = 0.0;
        double finalError = 0.0;
        std::size_t iterations = 0;
        bool converged = false;
    };

    BearingGaussNewtonResult OptimizeBearingTargetGaussNewton(const geodesy::GeoCoordinate &initialGuess, const std::vector<BearingObservation> &observations, const BearingGaussNewtonOptions &options = {});
}
