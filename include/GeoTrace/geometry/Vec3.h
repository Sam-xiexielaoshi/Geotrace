#pragma once

#include <cmath>

namespace geotrace::geometry
{
    struct Vec3
    {
        double x, y, z;
        double Length() const;
        Vec3 Normalize() const;
    };

    Vec3 Cross(const Vec3 &a, const Vec3 &b);
    double Dot(const Vec3 &a, const Vec3 &b);
}