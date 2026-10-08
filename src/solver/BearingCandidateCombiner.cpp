#include "GeoTrace/solver/BearingCandidateCombiner.h"

#include <stdexcept>

namespace geotrace::solver
{
    geometry::Vec3 CombineBearingCandidates(const std::vector<geometry::Vec3> &candidates)
    {
        if (candidates.empty())
        {
            throw std::invalid_argument("Cannot combine an empty set of bearing candidates.");
        }
        geometry::Vec3 combined{0.0, 0.0, 0.0};
        for (const auto &candidate : candidates)
        {
            combined = combined + candidate;
        }
        if (combined.IsNearlyZero())
        {
            throw std::runtime_error("Bearing candidates cancel out and do not define a target direction.");
        }
        return combined.Normalize();
    }
}