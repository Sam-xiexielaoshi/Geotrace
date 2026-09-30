#pragma once
#include "GeoTrace/geodesy/GeoCoordinate.h"
#include <string>
#include <vector>

namespace geotrace::io
{
    struct CsvResult
    {
        std::string caseId;
        geodesy::GeoCoordinate coordinate;
    };

    void WriteCsv(const std::string& filename, const std::vector<CsvResult>& results);
}