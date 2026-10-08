#include "GeoTrace/solver/BearingOptimizer.h"
#include "GeoTrace/solver/BearingObjective.h"

#include <algorithm>
#include <cmath>

namespace geotrace::solver
{
    geodesy::GeoCoordinate OptimizeBearingTarget(const geodesy::GeoCoordinate &initialGuess, const std::vector<BearingObservation> &observations)
    {
        double latitude = initialGuess.latitude;
        double longitude = initialGuess.longitude;
        double step = 1.0;
        while (step > 1e-5)
        {
            const geodesy::GeoCoordinate current{latitude, longitude};
            double bestError = BearingError(geodesy::LatLonToECEF(current), observations);
            bool improved = true;
            while (improved)
            {
                improved = false;

                const geodesy::GeoCoordinate candidates[] = {
                    {latitude + step, longitude},
                    {latitude - step, longitude},
                    {latitude, longitude + step},
                    {latitude, longitude - step}};

                for (const auto &candidate : candidates)
                {
                    if (candidate.latitude < -90.0 || candidate.latitude > 90.0)
                    {
                        continue;
                    }
                    double candidateLongitude = candidate.longitude;
                    while (candidateLongitude > 180.0)
                    {
                        candidateLongitude -= 360.0;
                    }
                    while (candidateLongitude < -180.0)
                    {
                        candidateLongitude += 360.0;
                    }

                    const geodesy::GeoCoordinate wrappedCandidate{candidate.latitude, candidateLongitude};
                    const double candidateError = BearingError(geodesy::LatLonToECEF(wrappedCandidate), observations);
                    if (candidateError < bestError)
                    {
                        latitude = wrappedCandidate.latitude;
                        longitude = wrappedCandidate.longitude;
                        bestError = candidateError;
                        improved = true;
                        break;
                    }
                }
            }
            step *= 0.5;
        }
        return {latitude, longitude};
    }
}