#include <catch2/catch_test_macros.hpp>

#include "GeoTrace/geometry/Vec3.h"

using namespace geotrace::geometry;

TEST_CASE("Vec3 length is calculated correctly")
{
    Vec3 v{3.0, 4.0, 0.0};

    REQUIRE(v.Length() == 5.0);
}

TEST_CASE("Vec3 normalization is calculated correctly")
{
    Vec3 v{3.0, 4.0, 0.0};

    Vec3 normalized = v.Normalize();

    REQUIRE(normalized.x == 0.6);
    REQUIRE(normalized.y == 0.8);
    REQUIRE(normalized.z == 0.0);
}

TEST_CASE("Vec3 cross product is calculated correctly")
{
    Vec3 xAxis{1.0, 0.0, 0.0};
    Vec3 yAxis{0.0, 1.0, 0.0};

    Vec3 result = Cross(xAxis, yAxis);

    REQUIRE(result.x == 0.0);
    REQUIRE(result.y == 0.0);
    REQUIRE(result.z == 1.0);
}

TEST_CASE("Vec3 dot product is calculated correctly")
{
    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{4.0, 5.0, 6.0};

    double result = Dot(a, b);

    REQUIRE(result == 32.0);
}

TEST_CASE("Vec3 addition is calculated correctly")
{
    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{4.0, 5.0, 6.0};

    Vec3 result = a + b;

    REQUIRE(result.x == 5.0);
    REQUIRE(result.y == 7.0);
    REQUIRE(result.z == 9.0);
}

TEST_CASE("Vec3 detects nearly zero vectors")
{
    const geotrace::geometry::Vec3 zero{0.0, 0.0, 0.0};

    REQUIRE(zero.IsNearlyZero());

    const geotrace::geometry::Vec3 tiny{
        1e-14,
        0.0,
        0.0};

    REQUIRE(tiny.IsNearlyZero());

    const geotrace::geometry::Vec3 normal{
        1.0,
        0.0,
        0.0};

    REQUIRE_FALSE(normal.IsNearlyZero());
}

TEST_CASE("Vec3 unary minus negates all components")
{
    const geotrace::geometry::Vec3 vector{
        1.0, -2.0, 3.0};

    const auto negated = -vector;

    REQUIRE(std::abs(negated.x + 1.0) < 1e-12);
    REQUIRE(std::abs(negated.y - 2.0) < 1e-12);
    REQUIRE(std::abs(negated.z + 3.0) < 1e-12);
}