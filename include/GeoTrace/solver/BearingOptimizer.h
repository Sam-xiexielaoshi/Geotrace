#pragma once

#include "GeoTrace/geodesy/GeoCoordinate.h"
#include "GeoTrace/solver/BearingObservation.h"

#include <vector>

namespace geotrace::solver
{
    geodesy::GeoCoordinate OptimizeBearingTarget(const geodesy::GeoCoordinate &initialGuess, const std::vector<BearingObservation> &observations);
}