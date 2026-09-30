#pragma once

#include "GeoTrace/geodesy/GeoCoordinate.h"

namespace geotrace::solver
{
    struct BearingObservation
    {
        geodesy::GeoCoordinate observer;
        double bearingDegrees;
    };
}