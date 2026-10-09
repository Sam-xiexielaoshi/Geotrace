
#include "GeoTrace/solver/BearingGaussNewtonOptimizer.h"

#include "GeoTrace/geodesy/GeoCoordinate.h"
#include "GeoTrace/solver/BearingJacobian.h"
#include "GeoTrace/solver/BearingObjective.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace geotrace::solver
{
    namespace
    {
        constexpr double kMinimumDamping = 1e-12;
        constexpr double kMaximumDamping = 1e12;
        constexpr int kMaximumLineSearchSteps = 12;

        double WrapLongitude(double longitude)
        {
            while (longitude > 180.0)
                longitude -= 360.0;
            while (longitude < -180.0)
                longitude += 360.0;
            return longitude;
        }

        double ComputeSquaredNorm(const std::vector<double> &values)
        {
            double result = 0.0;
            for (const double value : values)
                result += value * value;
            return result;
        }
    }

    BearingGaussNewtonResult OptimizeBearingTargetGaussNewton(const geodesy::GeoCoordinate &initialGuess, const std::vector<BearingObservation> &observations, const BearingGaussNewtonOptions &options)
    {
        if (observations.empty())
            throw std::invalid_argument("Gauss-Newton requires at least one observation.");

        if (!std::isfinite(initialGuess.latitude) || !std::isfinite(initialGuess.longitude) || initialGuess.latitude < -90.0 || initialGuess.latitude > 90.0)
        {
            throw std::invalid_argument("Initial guess must have valid finite coordinates.");
        }
        if (options.maxIterations == 0 || !(options.initialDamping > 0.0) ||
            !std::isfinite(options.initialDamping) || !(options.gradientTolerance > 0.0) ||
            !std::isfinite(options.gradientTolerance) || !(options.stepToleranceDegrees > 0.0) ||
            !std::isfinite(options.stepToleranceDegrees))
        {
            throw std::invalid_argument("Gauss-Newton options must be positive and finite.");
        }
        geodesy::GeoCoordinate current{initialGuess.latitude, WrapLongitude(initialGuess.longitude)};
        const auto objective = [&](const geodesy::GeoCoordinate &point)
        {
            return BearingError(geodesy::LatLonToECEF(point), observations);
        };
        BearingGaussNewtonResult result{};
        result.estimate = current;
        result.initialError = objective(current);
        result.finalError = result.initialError;
        double damping = options.initialDamping;
        for (std::size_t iteration = 0; iteration < options.maxIterations; ++iteration)
        {
            const std::vector<double> residuals = BearingResiduals(geodesy::LatLonToECEF(current), observations);
            const BearingJacobian jacobian = ComputeBearingJacobian(current, observations);
            // Build J^T J and J^T r.
            double a = 0.0;
            double b = 0.0;
            double c = 0.0;
            double g0 = 0.0;
            double g1 = 0.0;

            for (std::size_t row = 0; row < residuals.size(); ++row)
            {
                const double jLat = jacobian[row][0];
                const double jLon = jacobian[row][1];
                const double residual = residuals[row];
                a += jLat * jLat;
                b += jLat * jLon;
                c += jLon * jLon;
                g0 += jLat * residual;
                g1 += jLon * residual;
            }
            const double gradientNorm = std::hypot(g0, g1);
            if (!std::isfinite(gradientNorm))
                break;
            if (gradientNorm <= options.gradientTolerance)
            {
                result.converged = true;
                break;
            }
            bool accepted = false;
            for (int attempt = 0; attempt < kMaximumLineSearchSteps; ++attempt)
            {
                // Levenberg-style diagonal damping:
                // (J^T J + lambda D) delta = -J^T r.
                const double dampedA = a + damping * std::max(a, kMinimumDamping);
                const double dampedC = c + damping * std::max(c, kMinimumDamping);
                const double determinant = dampedA * dampedC - b * b;
                if (!(determinant > 0.0) || !std::isfinite(determinant))
                {
                    damping *= 10.0;
                    if (damping > kMaximumDamping)
                        break;
                    continue;
                }
                const double deltaLatitude = (-g0 * dampedC + b * g1) / determinant;
                const double deltaLongitude = (b * g0 - dampedA * g1) / determinant;
                if (!std::isfinite(deltaLatitude) || !std::isfinite(deltaLongitude))
                {
                    damping *= 10.0;
                    continue;
                }
                // Backtracking line search.
                double scale = 1.0;
                for (int lineSearch = 0; lineSearch < kMaximumLineSearchSteps; ++lineSearch)
                {
                    const double proposedLatitude = current.latitude + scale * deltaLatitude;
                    const double proposedLongitude = WrapLongitude(current.longitude + scale * deltaLongitude);
                    if (proposedLatitude >= -90.0 && proposedLatitude <= 90.0)
                    {
                        const geodesy::GeoCoordinate candidate{proposedLatitude, proposedLongitude};
                        const double candidateError = objective(candidate);
                        if (std::isfinite(candidateError) && candidateError < result.finalError)
                        {
                            current = candidate;
                            result.finalError = candidateError;
                            result.estimate = current;
                            result.iterations = iteration + 1;
                            accepted = true;
                            damping = std::max(damping * 0.3, kMinimumDamping);
                            const double stepSize = scale * std::hypot(deltaLatitude, deltaLongitude);
                            if (stepSize <= options.stepToleranceDegrees)
                            {
                                result.converged = true;
                            }
                            break;
                        }
                    }
                    scale *= 0.5;
                }
                if (accepted || result.converged)
                    break;
                damping *= 10.0;
                if (damping > kMaximumDamping)
                    break;
            }
            if (result.converged)
                break;
            // No acceptable step at the current damping.
            // Return the best estimate found so far.
            if (!accepted)
                break;
        }
        return result;
    }
}
