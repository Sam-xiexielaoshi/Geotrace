
#pragma once

#include "GeoTrace/solver/BearingObservation.h"
#include "GeoTrace/geodesy/GeoCoordinate.h"

#include <vector>

namespace geotrace::solver
{
    // Each row corresponds to an observation.
    // Column 0: derivative with respect to latitude (degrees).
    // Column 1: derivative with respect to longitude (degrees).
    using BearingJacobian = std::vector<std::vector<double>>;
    BearingJacobian ComputeBearingJacobian(const geodesy::GeoCoordinate &candidate, const std::vector<BearingObservation> &observations, double stepDegrees = 1e-3);
}
