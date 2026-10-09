
#include "GeoTrace/solver/BearingJacobian.h"

#include "GeoTrace/geodesy/GeoCoordinate.h"
#include "GeoTrace/solver/BearingObjective.h"

#include <stdexcept>

namespace geotrace::solver
{
    BearingJacobian ComputeBearingJacobian(const geodesy::GeoCoordinate &candidate, const std::vector<BearingObservation> &observations, double stepDegrees)
    {
        if (!(stepDegrees > 0.0))
        {
            throw std::invalid_argument("Jacobian step must be positive.");
        }
        const auto residualsAt = [&](const geodesy::GeoCoordinate &point)
        {
            return BearingResiduals(geodesy::LatLonToECEF(point), observations);
        };
        // const auto baseResiduals = residualsAt(candidate);
        BearingJacobian jacobian(observations.size(), std::vector<double>(2, 0.0));
        for (int column = 0; column < 2; ++column)
        {
            geodesy::GeoCoordinate minus = candidate;
            geodesy::GeoCoordinate plus = candidate;
            if (column == 0)
            {
                minus.latitude -= stepDegrees;
                plus.latitude += stepDegrees;
            }
            else
            {
                minus.longitude -= stepDegrees;
                plus.longitude += stepDegrees;
            }
            if (minus.latitude < -90.0 || plus.latitude > 90.0)
            {
                throw std::invalid_argument("Jacobian step crosses a latitude boundary.");
            }
            const auto minusResiduals = residualsAt(minus);
            const auto plusResiduals = residualsAt(plus);
            for (std::size_t row = 0; row < observations.size(); ++row)
            {
                jacobian[row][column] = (plusResiduals[row] - minusResiduals[row]) / (2.0 * stepDegrees);
            }
        }
        return jacobian;
    }
}
