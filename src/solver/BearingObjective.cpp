#include "GeoTrace/solver/BearingObjective.h"
#include "GeoTrace/geodesy/GeoCoordinate.h"

namespace geotrace::solver
{
    std::vector<double> BearingResiduals(const geometry::Vec3 &candidate, const std::vector<BearingObservation> &observations)
    {
        std::vector<double> residuals;
        residuals.reserve(observations.size());
        for (const auto &observation : observations)
        {
            const geometry::Vec3 observerPosition = geodesy::LatLonToECEF(observation.observer);
            const geometry::Vec3 bearingDirection = geodesy::BearingToDirection(observation.observer, observation.bearingDegrees);
            const geometry::Vec3 predictedTangent = geodesy::TargetTangentDirection(observerPosition, candidate);
            const double cosine = geometry::Dot(bearingDirection, predictedTangent);
            const double sine = geometry::Dot(observerPosition, geometry::Cross(bearingDirection, predictedTangent));
            residuals.push_back(std::atan2(sine, cosine));
        }
        return residuals;
    }

    double BearingError(const geometry::Vec3 &candidate, const std::vector<BearingObservation> &observations)
    {
        const std::vector<double> residuals = BearingResiduals(candidate, observations);
        double totalError = 0.0;
        for (const double residual : residuals)
        {
            totalError += residual * residual;
        }
        return totalError;
    }

}