#pragma once

#include "GeoTrace/geometry/BearingPlane.h"

namespace geotrace::geometry
{
    Vec3 IntersectPlanes(const BearingPlane& planeA, const BearingPlane& planeB);
}