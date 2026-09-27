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

    LocalFrame ComputeLocalFrame(const GeoCoordinate &coordinate)
    {
        const double latitude = DegreesToRadians(coordinate.latitude);
        const double longitude = DegreesToRadians(coordinate.longitude);
        const Vec3 north{
            -std::sin(latitude) * std::cos(longitude),
            -std::sin(latitude) * std::sin(longitude),
            std::cos(latitude)};

        const Vec3 east{
            -std::sin(longitude),
            std::cos(longitude),
            0.0};

        return {north, east};
    }

    Vec3 BearingToDirection(const GeoCoordinate &coordinate, double bearingDegrees)
    {
        const LocalFrame frame = ComputeLocalFrame(coordinate);
        const double bearning = DegreesToRadians(bearingDegrees);
        const double northWeight = std::cos(bearning);
        const double eastWeight = std::sin(bearning);
        return {
            frame.north.x * northWeight + frame.east.x * eastWeight,
            frame.north.y * northWeight + frame.east.y * eastWeight,
            frame.north.z * northWeight + frame.east.z * eastWeight};
    }
}