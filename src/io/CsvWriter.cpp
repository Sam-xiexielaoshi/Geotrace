#include "GeoTrace/io/CsvWriter.h"

#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace geotrace::io
{
    void WriteCsv(
        const std::string &filename,
        const std::vector<CsvResult> &results)
    {
        std::ofstream file(filename);

        if (!file.is_open())
        {
            throw std::runtime_error(
                "Failed to open file for writing: " + filename);
        }

        // Write CSV header
        file << "CaseID,Latitude,Longitude\n";

        // Write coordinates with six digits after the decimal point
        file << std::fixed << std::setprecision(6);

        for (const auto &result : results)
        {
            file << result.caseId << ","
                 << result.coordinate.latitude << ","
                 << result.coordinate.longitude
                 << "\n";
        }
    }
}