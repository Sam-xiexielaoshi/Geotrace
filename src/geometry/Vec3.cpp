#include "GeoTrace/geometry/Vec3.h"

namespace geotrace::geometry
{
    double Vec3::Length() const
    {
        return std::sqrt(x * x + y * y + z * z);
    }

    Vec3 Vec3::Normalize() const
    {
        const double length = Length();
        if (length == 0.0)
            return {0.0, 0.0, 0.0};
        return {x / length, y / length, z / length};
    }

    Vec3 Cross(const Vec3 &a, const Vec3 &b)
    {
        return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
    }

    double Dot(const Vec3 &a, const Vec3 &b)
    {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    Vec3 operator+(const Vec3 &a, const Vec3 &b)
    {
        return {a.x + b.x, a.y + b.y, a.z + b.z};
    }
}