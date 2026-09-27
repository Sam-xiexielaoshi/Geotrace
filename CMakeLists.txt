cmake_minimum_required(VERSION 3.25)

project(GeoTrace
    VERSION 0.1.0
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(GeoTrace
    src/main.cpp
    src/geometry/Vec3.cpp
    src/geodesy/GeoCoordinate.cpp
    src/geometry/BearingPlane.cpp
    src/geometry/PlaneIntersection.cpp
)

target_include_directories(GeoTrace
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)