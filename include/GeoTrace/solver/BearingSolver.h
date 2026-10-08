#pragma once

#include "GeoTrace/solver/BearingObservation.h"
#include "GeoTrace/geodesy/GeoCoordinate.h"

#include <vector>

namespace geotrace::solver
{
    geodesy::GeoCoordinate SolveInitial(const std::vector<BearingObservation>& observations);
    geodesy::GeoCoordinate Solve(const std::vector<BearingObservation>& observations);
}