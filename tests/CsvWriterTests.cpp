#include <catch2/catch_test_macros.hpp>

#include "GeoTrace/io/CsvWriter.h"

#include <fstream>
#include <sstream>
#include <string>
#include <cstdio>

TEST_CASE("CSV writer writes a single result correctly")
{
    const std::string filename = "test_output.csv";

    const geotrace::io::CsvResult result{
        "CASE001",
        {20.0, 40.0}};

    geotrace::io::WriteCsv(
        filename,
        {result});

    std::ifstream file(filename);

    REQUIRE(file.is_open());

    std::stringstream contents;
    contents << file.rdbuf();

    const std::string output = contents.str();

    REQUIRE(
        output ==
        "CaseID,Latitude,Longitude\n"
        "CASE001,20.000000,40.000000\n");

    file.close();

    std::remove(filename.c_str());
}