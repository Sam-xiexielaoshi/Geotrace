#include "GeoTrace/geodesy/GeoCoordinate.h"
#include "GeoTrace/geometry/Vec3.h"

#include <cmath>

namespace geotrace::geodesy
{
    using geotrace::geometry::Vec3;
    constexpr double PI = 3.14159265358979323846;
    double DegreesToRadians(double degrees)
    {
        return degrees * PI / 180.0;
    }

    Vec3 LatLonToECEF(const GeoCoordinate &coordinate)
    {
        const double latitude = DegreesToRadians(coordinate.latitude);

        const double longitude = DegreesToRadians(coordinate.longitude);

        return {
            std::cos(latitude) * std::cos(longitude),
            std::cos(latitude) * std::sin(longitude),
            std::sin(latitude)};
    }
}