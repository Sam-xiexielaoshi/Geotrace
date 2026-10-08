#pragma once

#include "GeoTrace/geometry/Vec3.h"

#include <vector>

namespace geotrace::solver
{
    geometry::Vec3 CombineBearingCandidates(const std::vector<geometry::Vec3> &candidates);
}