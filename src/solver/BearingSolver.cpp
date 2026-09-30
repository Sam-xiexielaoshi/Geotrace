#include "GeoTrace/solver/BearingSolver.h"

#include "GeoTrace/geometry/BearingPlane.h"
#include "GeoTrace/geometry/PlaneIntersection.h"
#include "GeoTrace/geodesy/GeoCoordinate.h"

namespace geotrace::solver
{
    geodesy::GeoCoordinate Solve(const BearingObservation &observationA, const BearingObservation &observationB, const BearingObservation &observationC)
    {
        // convert observer positions to 3d unit vectors
        const auto positionA = geodesy::LatLonToECEF(observationA.observer);
        const auto positionB = geodesy::LatLonToECEF(observationB.observer);
        const auto positionC = geodesy::LatLonToECEF(observationC.observer);

        // convert bearing tino 3d tangent directions
        const auto directionA = geodesy::BearingToDirection(observationA.observer, observationA.bearingDegrees);
        const auto directionB = geodesy::BearingToDirection(observationB.observer, observationB.bearingDegrees);
        const auto directionC = geodesy::BearingToDirection(observationC.observer, observationC.bearingDegrees);

        // construct bearing plane
        const auto planeA = geometry::CreateBearingPlane(positionA, directionA);
        const auto planeB = geometry::CreateBearingPlane(positionB, directionB);
        const auto planeC = geometry::CreateBearingPlane(positionC, directionC);

        // calculate pairwise plane intersections
        const auto intersectionAB = geometry::IntersectPlanes(planeA, planeB);
        const auto intersectionBC = geometry::IntersectPlanes(planeB, planeC);
        const auto intersectionCA = geometry::IntersectPlanes(planeA, planeC);

        // combine the intersections
        const auto combined = (intersectionAB + intersectionBC + intersectionCA).Normalize();

        //convert final direction back to latitude and longitude
        return geodesy::ECEFToLatLon(combined);
    }
}