#pragma once

#include "GeoTrace/geometry/BearingPlane.h"
#include "GeoTrace/solver/BearingObservation.h"

namespace geotrace::geometry
{
    Vec3 IntersectPlanes(const BearingPlane &planeA, const BearingPlane &planeB);
    Vec3 ResolveIntersectionDirection(const Vec3 &intersection, const solver::BearingObservation &observationA, const solver::BearingObservation &observationB);
}