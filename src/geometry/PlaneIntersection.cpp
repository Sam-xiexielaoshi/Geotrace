#include "GeoTrace/geometry/PlaneIntersection.h"

namespace geotrace::geometry
{
    Vec3 IntersectPlanes(const BearingPlane& planeA, const BearingPlane& planeB)
    {
        return Cross(planeA.normal, planeB.normal).Normalize();
    }
}