#pragma once

#include "GeoTrace/geometry/Vec3.h"

namespace geotrace::geometry
{
    struct BearingPlane
    {
        Vec3 normal;
    };

    BearingPlane CreateBearingPlane(const Vec3& observerPosition, const Vec3& bearingDirection);
}