#pragma once
#include "GeoTrace/io/CsvRecord.h"
#include <string>
#include <vector>

namespace geotrace::io
{
    std::vector<CsvRecord> ReadCsv(const std::string& filename);
}