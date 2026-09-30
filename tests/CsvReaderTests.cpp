#include <catch2/catch_test_macros.hpp>
#include "GeoTrace/io/CsvReader.h"

#include <fstream>
#include <string>

TEST_CASE("CSV reader parses a three-observer record")
{
    const std::string filename = "test_input.csv";

    {
        std::ofstream file(filename);

        file << "CaseID,"
             << "ObserverA_Lat,ObserverA_Lon,ObserverA_Bearing,"
             << "ObserverB_Lat,ObserverB_Lon,ObserverB_Bearing,"
             << "ObserverC_Lat,ObserverC_Lon,ObserverC_Bearing\n";

        file << "CASE001,"
             << "0.0,0.0,60.479848,"
             << "0.0,90.0,295.413767,"
             << "30.0,45.0,205.480025\n";
    }

    const auto records = geotrace::io::ReadCsv(filename);

    REQUIRE(records.size() == 1);

    const auto &record = records[0];

    REQUIRE(record.caseId == "CASE001");

    REQUIRE(record.observationA.observer.latitude == 0.0);
    REQUIRE(record.observationA.observer.longitude == 0.0);
    REQUIRE(record.observationA.bearingDegrees == 60.479848);

    REQUIRE(record.observationB.observer.latitude == 0.0);
    REQUIRE(record.observationB.observer.longitude == 90.0);
    REQUIRE(record.observationB.bearingDegrees == 295.413767);

    REQUIRE(record.observationC.observer.latitude == 30.0);
    REQUIRE(record.observationC.observer.longitude == 45.0);
    REQUIRE(record.observationC.bearingDegrees == 205.480025);

    std::remove(filename.c_str());
}

TEST_CASE("CSV reader throws when file does not exist")
{
    const std::string filename = "this_file_does_not_exist.csv";

    REQUIRE_THROWS(
        geotrace::io::ReadCsv(filename)
    );
}

TEST_CASE("CSV reader throws when a row has the wrong number of fields")
{
    const std::string filename = "invalid_column_count.csv";

    {
        std::ofstream file(filename);

        file << "CaseID,"
             << "ObserverA_Lat,ObserverA_Lon,ObserverA_Bearing,"
             << "ObserverB_Lat,ObserverB_Lon,ObserverB_Bearing,"
             << "ObserverC_Lat,ObserverC_Lon,ObserverC_Bearing\n";

        // Only 9 fields instead of the required 10.
        file << "CASE001,"
             << "0.0,0.0,60.479848,"
             << "0.0,90.0,295.413767,"
             << "30.0,45.0\n";
    }

    REQUIRE_THROWS(
        geotrace::io::ReadCsv(filename)
    );

    std::remove(filename.c_str());
}

TEST_CASE("CSV reader throws when a numeric field is invalid")
{
    const std::string filename = "invalid_numeric_value.csv";

    {
        std::ofstream file(filename);

        file << "CaseID,"
             << "ObserverA_Lat,ObserverA_Lon,ObserverA_Bearing,"
             << "ObserverB_Lat,ObserverB_Lon,ObserverB_Bearing,"
             << "ObserverC_Lat,ObserverC_Lon,ObserverC_Bearing\n";

        // "INVALID" cannot be converted to a double.
        file << "CASE001,"
             << "INVALID,0.0,60.479848,"
             << "0.0,90.0,295.413767,"
             << "30.0,45.0,205.480025\n";
    }

    REQUIRE_THROWS(
        geotrace::io::ReadCsv(filename)
    );

    std::remove(filename.c_str());
}