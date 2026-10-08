#pragma once

#include "GeoTrace/geometry/BearingPlane.h"
#include "GeoTrace/solver/BearingObservation.h"

#include <vector>

namespace geotrace::solver
{
    std::vector<geometry::BearingPlane> BuildBearingPlanes(const std::vector<BearingObservation> &observations);
}