#include "GeoTrace/geometry/PlaneIntersection.h"
#include "GeoTrace/geodesy/GeoCoordinate.h"
#include <stdexcept>

namespace geotrace::geometry
{
    Vec3 IntersectPlanes(const BearingPlane &planeA, const BearingPlane &planeB)
    {
        const Vec3 intersection = Cross(planeA.normal, planeB.normal);

        if (intersection.IsNearlyZero())
        {
            throw std::runtime_error("Planes are parallel and do not intersect.");
        }
        return intersection.Normalize();
    }

    Vec3 ResolveIntersectionDirection(const Vec3 &intersection, const solver::BearingObservation &observationA, const solver::BearingObservation &observationB)
    {
        const Vec3 candidateA = intersection;
        const Vec3 candidateB = -intersection;

        const Vec3 positionA = geodesy::LatLonToECEF(observationA.observer);
        const Vec3 bearingDirectionA = geodesy::BearingToDirection(observationA.observer, observationA.bearingDegrees);

        const Vec3 positionB = geodesy::LatLonToECEF(observationB.observer);
        const Vec3 bearingDirectionB = geodesy::BearingToDirection(observationB.observer, observationB.bearingDegrees);

        const bool candidateAConsistent = geodesy::IsBearingConsistent(positionA, bearingDirectionA, candidateA) && geodesy::IsBearingConsistent(positionB, bearingDirectionB, candidateA);
        const bool candidateBConsistent = geodesy::IsBearingConsistent(positionA, bearingDirectionA, candidateB) && geodesy::IsBearingConsistent(positionB, bearingDirectionB, candidateB);

        if (candidateAConsistent && !candidateBConsistent)
        {
            return candidateA;
        }
        if (candidateBConsistent && !candidateAConsistent)
        {
            return candidateB;
        }
        throw std::runtime_error("Could not resolve antipodal intersection direction");
    }
}