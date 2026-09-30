#pragma once

#include "GeoTrace/solver/BearingObservation.h"
#include "GeoTrace/geodesy/GeoCoordinate.h"

namespace geotrace::solver
{
    geodesy::GeoCoordinate Solve(const BearingObservation &observationA, const BearingObservation &observationB, const BearingObservation &observationC);
}