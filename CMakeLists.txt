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
    src/solver/BearingSolver.cpp
    src/io/CsvReader.cpp
    src/io/CsvWriter.cpp
    src/solver/BearingObjective.cpp
    src/solver/BearingOptimizer.cpp
    src/solver/BearingPlaneBuilder.cpp
    src/solver/BearingIntersectionBuilder.cpp
    src/solver/BearingCandidateCombiner.cpp
    src/solver/BearingJacobian.cpp
    src/solver/BearingGaussNewtonOptimizer.cpp

)

target_include_directories(GeoTrace

    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include

)

include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.8.1
)

FetchContent_MakeAvailable(Catch2)

enable_testing()

add_executable(GeoTraceTests

    tests/Vec3Tests.cpp
    tests/GeoCoordinateTests.cpp
    tests/BearingPlaneTests.cpp
    tests/PlaneIntersectionTests.cpp
    tests/BearingSolverTests.cpp
    tests/CsvReaderTests.cpp
    tests/CsvWriterTests.cpp

    src/geometry/Vec3.cpp
    src/geodesy/GeoCoordinate.cpp
    src/geometry/BearingPlane.cpp
    src/geometry/PlaneIntersection.cpp
    src/solver/BearingSolver.cpp
    src/io/CsvReader.cpp
    src/io/CsvWriter.cpp
    src/solver/BearingObjective.cpp
    src/solver/BearingOptimizer.cpp
    src/solver/BearingPlaneBuilder.cpp
    src/solver/BearingIntersectionBuilder.cpp
    src/solver/BearingCandidateCombiner.cpp
    src/solver/BearingJacobian.cpp
    src/solver/BearingGaussNewtonOptimizer.cpp

)

target_include_directories(GeoTraceTests

    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include

)

target_link_libraries(GeoTraceTests

    PRIVATE
        Catch2::Catch2WithMain

)

include(CTest)
include(Catch)

catch_discover_tests(GeoTraceTests)