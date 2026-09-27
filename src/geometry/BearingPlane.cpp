#include "GeoTrace/geometry/BearingPlane.h"

namespace geotrace::geometry
{
    BearingPlane CreateBearingPlane(const Vec3& observerPosition, const Vec3& bearingDirection)
    {
        return {Cross(observerPosition, bearingDirection).Normalize()};
    }
}