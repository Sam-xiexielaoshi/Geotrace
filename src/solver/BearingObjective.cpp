#include "GeoTrace/solver/BearingObjective.h"
#include "GeoTrace/geodesy/GeoCoordinate.h"

namespace geotrace::solver
{
    double BearingError(const geometry::Vec3& candidate, const std::vector<BearingObservation>& observations)
    {
        double totalError = 0.0;
        for (const auto& observation : observations)
        {
            const geometry::Vec3 observerPosition = geodesy::LatLonToECEF(observation.observer);
            const geometry::Vec3 bearingDirection = geodesy::BearingToDirection(observation.observer, observation.bearingDegrees);
            const double residual = geodesy::BearingAngularResidual(observerPosition, bearingDirection, candidate);
            totalError += residual * residual;
        }
        return totalError;
    }
}