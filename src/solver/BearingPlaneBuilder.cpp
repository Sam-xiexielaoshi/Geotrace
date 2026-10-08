#include "GeoTrace/solver/BearingPlaneBuilder.h"
#include "GeoTrace/geodesy/GeoCoordinate.h"
#include <stdexcept>

namespace geotrace::solver
{
    std::vector<geometry::BearingPlane> BuildBearingPlanes(const std::vector<BearingObservation> &observations)
    {
        if (observations.size() < 3)
        {
            throw std::invalid_argument("At least 3 bearing observations are required.");
        }
        std::vector<geometry::BearingPlane> planes;
        planes.reserve(observations.size());
        for (const auto &observation : observations)
        {
            const auto position = geodesy::LatLonToECEF(observation.observer);
            const auto direction = geodesy::BearingToDirection(observation.observer, observation.bearingDegrees);
            const auto plane = geometry::CreateBearingPlane(position, direction);
            planes.push_back(plane);
        }
        return planes;
    }
}