#pragma once

#include "GeoTrace/geometry/Vec3.h"
#include "GeoTrace/geometry/BearingPlane.h"
#include "GeoTrace/solver/BearingObservation.h"

#include <vector>

namespace geotrace::solver
{
    std::vector<geometry::Vec3> BuildBearingIntersections(const std::vector<geometry::BearingPlane> &planes, const std::vector<BearingObservation> &observations);
}