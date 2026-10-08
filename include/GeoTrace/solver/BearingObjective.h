#pragma once
#include "GeoTrace/geometry/Vec3.h"
#include "GeoTrace/solver/BearingObservation.h"

#include <vector>

namespace geotrace::solver
{
    double BearingError(const geometry::Vec3 &candidate, const std::vector<BearingObservation> &observations);
}