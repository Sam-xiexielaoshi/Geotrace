#include <iostream>
#include <iomanip>
#include "GeoTrace/geometry/Vec3.h"
#include "GeoTrace/geodesy/GeoCoordinate.h"
#include "GeoTrace/geometry/BearingPlane.h"
#include "GeoTrace/geometry/PlaneIntersection.h"

using namespace geotrace::geometry;
using namespace geotrace::geodesy;

void PrintECEF(const char *name, const GeoCoordinate &coordinate)
{
    const auto position = LatLonToECEF(coordinate);
    std::cout << name << "\n";
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

void PrintBearingDirection(const char *name, const GeoCoordinate &coordinate, double bearingDegrees)
{
    const auto direction = BearingToDirection(coordinate, bearingDegrees);
    std::cout << name << " | Bearing: " << bearingDegrees << " \n";
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
    std::cout << "\nLength of a: " << a.Length() << "\n";

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

    const auto position = LatLonToECEF({0.0, 0.0});
    const auto direction = BearingToDirection({0.0, 0.0}, 45.0);
    const auto plane = CreateBearingPlane(position, direction);

    std::cout << "Bearing plane normal: (" << plane.normal.x << ", " << plane.normal.y << ", " << plane.normal.z << ")\n";
    std::cout << "|Normal| = " << plane.normal.Length() << "\n";
    std::cout << "Position . Normal = " << Dot(position, plane.normal) << "\n";
    std::cout << "Direction . Normal = " << Dot(direction, plane.normal) << "\n";

    BearingPlane planeA{
        {1.0, 0.0, 0.0}};

    BearingPlane planeB{
        {0.0, 1.0, 0.0}};

    Vec3 intersection = IntersectPlanes(planeA, planeB);
    std::cout << "Plane intersection: (" << intersection.x << ", " << intersection.y << ", " << intersection.z << ")\n";
    std::cout << "|Intersection| = " << intersection.Length() << "\n";
    std::cout << "Intersection . PlaneA = " << Dot(intersection, planeA.normal) << "\n";
    std::cout << "Intersection . PlaneB = " << Dot(intersection, planeB.normal) << "\n";

    // ------------------------------------------------------------
    // End-to-end reconstruction test
    // ------------------------------------------------------------

    const GeoCoordinate target{20.0, 40.0};

    const GeoCoordinate observerA{0.0, 0.0};
    const GeoCoordinate observerB{0.0, 90.0};

    // Generate synthetic bearings from the known target.
    // In the real application these bearings would come from input data.
    const double bearingA = InitialBearing(observerA, target);
    const double bearingB = InitialBearing(observerB, target);

    std::cout << "\n";
    std::cout << "========== End-to-End Reconstruction ==========\n";

    std::cout << "Known target: ("
              << target.latitude << ", "
              << target.longitude << ")\n";

    std::cout << "Observer A bearing: "
              << bearingA << " degrees\n";

    std::cout << "Observer B bearing: "
              << bearingB << " degrees\n";

    // ------------------------------------------------------------
    // Convert observers to 3D positions
    // ------------------------------------------------------------

    const Vec3 positionA = LatLonToECEF(observerA);
    const Vec3 positionB = LatLonToECEF(observerB);

    // ------------------------------------------------------------
    // Convert bearings into 3D tangent directions
    // ------------------------------------------------------------

    const Vec3 directionA = BearingToDirection(observerA, bearingA);
    const Vec3 directionB = BearingToDirection(observerB, bearingB);

    // ------------------------------------------------------------
    // Construct the two bearing planes
    // ------------------------------------------------------------

    const BearingPlane bearingPlaneA =
        CreateBearingPlane(positionA, directionA);

    const BearingPlane bearingPlaneB =
        CreateBearingPlane(positionB, directionB);

    // ------------------------------------------------------------
    // Intersect the two bearing planes
    // ------------------------------------------------------------

    const Vec3 targetIntersection =
        IntersectPlanes(bearingPlaneA, bearingPlaneB);

    // The intersection is a line through the Earth's center,
    // so both directions are possible.
    const Vec3 candidateA = targetIntersection;
    const Vec3 candidateB = {
        -targetIntersection.x,
        -targetIntersection.y,
        -targetIntersection.z};

    // ------------------------------------------------------------
    // Convert both candidates back to latitude / longitude
    // ------------------------------------------------------------

    const GeoCoordinate resultA =
        ECEFToLatLon(candidateA);

    const GeoCoordinate resultB =
        ECEFToLatLon(candidateB);

    // ------------------------------------------------------------
    // Print results
    // ------------------------------------------------------------

    std::cout << "\nCandidate 1: ("
              << resultA.latitude << ", "
              << resultA.longitude << ")\n";

    std::cout << "Candidate 2: ("
              << resultB.latitude << ", "
              << resultB.longitude << ")\n";

    std::cout << "===============================================\n";
    return 0;
}