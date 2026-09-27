#pragma once
#include "GeoTrace/geometry/Vec3.h"

namespace geotrace::geodesy
{
    struct GeoCoordinate
    {
        double latitude;
        double longitude;
    };

    struct LocalFrame
    {
        geometry::Vec3 north;
        geometry::Vec3 east;
    };

    double DegreesToRadians(double degrees);

    geometry::Vec3 LatLonToECEF(const GeoCoordinate &coordinate);
    LocalFrame ComputeLocalFrame(const GeoCoordinate &coordinate);
    geometry::Vec3 BearingToDirection(const GeoCoordinate &coordinate, double bearingDegrees);
}