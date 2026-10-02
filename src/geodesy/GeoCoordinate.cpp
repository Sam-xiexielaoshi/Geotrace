#include "GeoTrace/geodesy/GeoCoordinate.h"
#include "GeoTrace/geometry/Vec3.h"

#include <cmath>
#include <stdexcept>

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
        const double bearing = DegreesToRadians(bearingDegrees);
        const double northWeight = std::cos(bearing);
        const double eastWeight = std::sin(bearing);
        return {
            frame.north.x * northWeight + frame.east.x * eastWeight,
            frame.north.y * northWeight + frame.east.y * eastWeight,
            frame.north.z * northWeight + frame.east.z * eastWeight};
    }

    double InitialBearing(const GeoCoordinate &observer, const GeoCoordinate &target)
    {
        const double latitude1 = DegreesToRadians(observer.latitude);
        const double latitude2 = DegreesToRadians(target.latitude);
        const double longitude1 = DegreesToRadians(observer.longitude);
        const double longitude2 = DegreesToRadians(target.longitude);
        const double deltaLongitude = longitude2 - longitude1;
        const double y = std::sin(deltaLongitude) * std::cos(latitude2);
        const double x = std::cos(latitude1) * std::sin(latitude2) - std::sin(latitude1) * std::cos(latitude2) * std::cos(deltaLongitude);
        double bearing = std::atan2(y, x);
        bearing = bearing * 180.0 / PI;
        if (bearing < 0.0)
            bearing += 360.0;
        return bearing;
    }

    GeoCoordinate ECEFToLatLon(const geometry::Vec3 &position)
    {
        const double latitude = std::asin(position.z);
        const double longitude = std::atan2(position.y, position.x);
        return {latitude * 180.0 / PI, longitude * 180.0 / PI};
    }

    geometry::Vec3 TargetTangentDirection(const geometry::Vec3 &observer, const geometry::Vec3 &target)
    {
        const double projection = geometry::Dot(observer, target);
        const geometry::Vec3 tangent{target.x - observer.x * projection, target.y - observer.y * projection, target.z - observer.z * projection};

        if (tangent.IsNearlyZero())
        {
            throw std::runtime_error("Observer and target are in the same direction; tangent direction is undefined.");
        }
        return target.Normalize();
    }
}