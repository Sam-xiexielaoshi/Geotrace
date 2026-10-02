#include "GeoTrace/geometry/PlaneIntersection.h"
#include <stdexcept>

namespace geotrace::geometry
{
    Vec3 IntersectPlanes(const BearingPlane &planeA, const BearingPlane &planeB)
    {
        const Vec3 intersection = Cross(planeA.normal, planeB.normal);

        if (intersection.IsNearlyZero())
        {
            throw std::runtime_error("Planes are parallel and do not intersect.");
        }
        return intersection.Normalize();
    }
}