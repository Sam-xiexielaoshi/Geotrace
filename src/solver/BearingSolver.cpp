#include "GeoTrace/solver/BearingSolver.h"
#include "GeoTrace/solver/BearingOptimizer.h"

#include "GeoTrace/solver/BearingPlaneBuilder.h"
#include "GeoTrace/solver/BearingIntersectionBuilder.h"
#include "GeoTrace/solver/BearingCandidateCombiner.h"

#include <stdexcept>
#include <vector>

namespace geotrace::solver
{
    geodesy::GeoCoordinate SolveInitial(const std::vector<BearingObservation> &observations)
    {
        if (observations.size() < 3)
        {
            throw std::invalid_argument(
                "At least 3 bearing observations are required.");
        }
        const auto planes = BuildBearingPlanes(observations);
        const auto intersections = BuildBearingIntersections(planes, observations);
        const auto combined = CombineBearingCandidates(intersections);
        return geodesy::ECEFToLatLon(combined);
    }

    geodesy::GeoCoordinate Solve(const std::vector<BearingObservation> &observations)
    {
        const auto initialEstimate = SolveInitial(observations);
        return OptimizeBearingTarget(initialEstimate, observations);
    }
}