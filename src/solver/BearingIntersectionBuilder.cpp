#include "GeoTrace/solver/BearingIntersectionBuilder.h"
#include "GeoTrace/geometry/PlaneIntersection.h"
#include <stdexcept>

namespace geotrace::solver
{
    std::vector<geometry::Vec3> BuildBearingIntersections(const std::vector<geometry::BearingPlane> &planes, const std::vector<BearingObservation> &observations)
    {
        if (planes.size() != observations.size())
        {
            throw std::invalid_argument("The number of planes must match the number of observations.");
        }
        if (planes.size() < 3)
        {
            throw std::invalid_argument(
                "At least 3 bearing planes are required.");
        }
        std::vector<geometry::Vec3> intersections;
        const std::size_t observerCount = planes.size();
        const std::size_t pairCount = observerCount * (observerCount - 1) / 2;
        intersections.reserve(pairCount);

        for (size_t i = 0; i < observerCount; ++i)
        {
            for (size_t j = i + 1; j < observerCount; ++j)
            {
                const auto rawIntersection = geometry::IntersectPlanes(planes[i], planes[j]);
                const auto resolvedIntersection = geometry::ResolveIntersectionDirection(rawIntersection, observations[i], observations[j]);
                intersections.push_back(resolvedIntersection);
            }
        }

        return intersections;
    }
}