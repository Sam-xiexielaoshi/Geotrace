#include "GeoTrace/solver/BearingSolver.h"
#include "GeoTrace/geometry/BearingPlane.h"
#include "GeoTrace/solver/BearingOptimizer.h"

#include "GeoTrace/geometry/PlaneIntersection.h"
#include "GeoTrace/geodesy/GeoCoordinate.h"

#include <stdexcept>
#include <vector>

namespace geotrace::solver
{
    geodesy::GeoCoordinate SolveInitial(const BearingObservation &observationA, const BearingObservation &observationB, const BearingObservation &observationC)
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
        const auto rawIntersectionAB = geometry::IntersectPlanes(planeA, planeB);
        const auto rawIntersectionBC = geometry::IntersectPlanes(planeB, planeC);
        const auto rawIntersectionCA = geometry::IntersectPlanes(planeC, planeA);
        const auto intersectionAB = geometry::ResolveIntersectionDirection(rawIntersectionAB, observationA, observationB);
        const auto intersectionBC = geometry::ResolveIntersectionDirection(rawIntersectionBC, observationB, observationC);
        const auto intersectionCA = geometry::ResolveIntersectionDirection(rawIntersectionCA, observationC, observationA);

        // combine the intersections
        const auto combined = intersectionAB + intersectionBC + intersectionCA;
        if (combined.IsNearlyZero())
        {
            throw std::runtime_error("Bearing observations do not converge to a single point");
        }

        // convert final direction back to latitude and longitude
        return geodesy::ECEFToLatLon(combined.Normalize());
    }

    geodesy::GeoCoordinate Solve(const BearingObservation &observationA, const BearingObservation &observationB, const BearingObservation &observationC)
    {
        const auto initialEstimate = SolveInitial(observationA, observationB, observationC);
        const std::vector<BearingObservation> observations{observationA, observationB, observationC};
        return OptimizeBearingTarget(initialEstimate, observations);
    }
}