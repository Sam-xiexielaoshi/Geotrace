#pragma once
#include "GeoTrace/geometry/Vec3.h"

namespace geotrace::geodesy
{
    struct GeoCoordinate
    {
        double latitude;
        double longitude;
    };

    double DegreesToRadians(double degrees);

    geometry::Vec3 LatLonToECEF(const GeoCoordinate &coordinate);
}