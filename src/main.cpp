#include <iostream>
#include <iomanip>
#include "GeoTrace/geometry/Vec3.h"
#include "GeoTrace/geodesy/GeoCoordinate.h"

using namespace geotrace::geometry;
using namespace geotrace::geodesy;

void PrintECEF(const char *name, const GeoCoordinate &coordinate)
{
    const auto position = LatLonToECEF(coordinate);
    std::cout << name << "/n";
    std::cout << "ECEF: (" << position.x << ", " << position.y << ", " << position.z << ")\n";
    std::cout << "\n";
}

void PrintLocalFrame(const char *name, const GeoCoordinate &coordinate)
{
    const auto frame = ComputeLocalFrame(coordinate);
    std::cout << name << '\n';
    std::cout << "North: " << frame.north.x << ", " << frame.north.y << ", " << frame.north.z << "\n";
    std::cout << "East: " << frame.east.x << ", " << frame.east.y << ", " << frame.east.z << "\n";

    std::cout
        << "|North| = "
        << frame.north.Length()
        << '\n';

    std::cout
        << "|East| = "
        << frame.east.Length()
        << '\n';
}

void PrintBearingDirection(const char *name, const GeoCoordinate &coordinate, double bearningDegrees)
{
    const auto direction = BearingToDirection(coordinate, bearningDegrees);
    std::cout << name << " | Bearing: " << bearningDegrees << " \n";
    std::cout << "Direction: " << direction.x << ", " << direction.y << ", " << direction.z << "\n";
    std::cout << "|Direction| = " << direction.Length() << "\n";
}

int main()
{
    std::cout << std::fixed << std::setprecision(6);

    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{4.0, 5.0, 6.0};

    Vec3 crossProduct = Cross(a, b);
    std::cout << "Cross product: (" << crossProduct.x << ", " << crossProduct.y << ", " << crossProduct.z << ")\n";
    std::cout << "\nLenght of a: " << a.Length() << "\n";

    PrintECEF("Equator / Prime Meridian", {0.0, 0.0});
    PrintECEF("Equator / 90E", {0.0, 90.0});
    PrintECEF("North Pole", {90.0, 0.0});

    PrintLocalFrame("Equator / Prime Meridian", {0.0, 0.0});
    PrintLocalFrame("Equator / 90E", {0.0, 90.0});
    PrintLocalFrame("North Pole", {90.0, 0.0});
    PrintBearingDirection("Equator / Prime Meridian", {0.0, 0.0}, 0.0);
    PrintBearingDirection("Equator / Prime Meridian", {0.0, 0.0}, 90.0);
    PrintBearingDirection("Equator / Prime Meridian", {0.0, 0.0}, 180.0);
    PrintBearingDirection("Equator / Prime Meridian", {0.0, 0.0}, 270.0);
    PrintBearingDirection("Equator / Prime Meridian", {0.0, 0.0}, 45.0);
    return 0;
}