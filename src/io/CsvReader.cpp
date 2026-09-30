#include "GeoTrace/io/CsvReader.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace geotrace::io
{
    namespace
    {
        std::vector<std::string> SplitCsvLine(const std::string &line)
        {
            std::vector<std::string> fields;
            std::stringstream stream(line);
            std::string field;
            while (std::getline(stream, field, ','))
            {
                fields.push_back(field);
            }
            return fields;
        }
        double ParseDouble(const std::string &value)
        {
            try
            {
                return std::stod(value);
            }
            catch (...)
            {
                throw ::std::runtime_error("Failed to parse double from string: " + value);
            }
        }
    }

    std::vector<CsvRecord> ReadCsv(const std::string &filename)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {
            throw std::runtime_error("Failed to open CSV file: " + filename);
        }
        std::vector<CsvRecord> records;
        std::string line;
        // spiking the header
        std::getline(file, line);
        while (std::getline(file, line))
        {
            if (line.empty())
            {
                continue; // Skip empty lines
            }
            const auto fields = SplitCsvLine(line);
            if (fields.size() != 10)
            {
                throw std::runtime_error("Expected 10 CSV fields, got " +
                                         std::to_string(fields.size()));
            }
            CsvRecord record;
            record.caseId = fields[0];
            record.observationA.observer = {ParseDouble(fields[1]), ParseDouble(fields[2])};
            record.observationA.bearingDegrees = ParseDouble(fields[3]);
            record.observationB.observer = {ParseDouble(fields[4]), ParseDouble(fields[5])};
            record.observationB.bearingDegrees = ParseDouble(fields[6]);
            record.observationC.observer = {ParseDouble(fields[7]), ParseDouble(fields[8])};
            record.observationC.bearingDegrees = ParseDouble(fields[9]);

            records.push_back(record);
        }
        return records;
    }
}
