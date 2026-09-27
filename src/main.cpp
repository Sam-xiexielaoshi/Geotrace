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

int main()
{
    std::cout<<std::fixed<<std::setprecision(6);

    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{4.0, 5.0, 6.0};

    Vec3 crossProduct = Cross(a, b);
    std::cout << "Cross product: (" << crossProduct.x << ", " << crossProduct.y << ", " << crossProduct.z << ")\n";
    std::cout << "\nLenght of a: " << a.Length() << "\n";

    PrintECEF("Equator / Prime Meridian", {0.0, 0.0});
    PrintECEF("Equator / 90E", {0.0, 90.0});
    PrintECEF("North Pole", {90.0, 0.0});
    return 0;
}