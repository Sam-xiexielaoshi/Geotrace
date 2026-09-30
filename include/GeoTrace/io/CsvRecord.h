#pragma once
#include "GeoTrace/solver/BearingObservation.h"
#include <string>

namespace geotrace::io
{
    struct CsvRecord
    {
        std::string caseId;

        solver::BearingObservation observationA;
        solver::BearingObservation observationB;
        solver::BearingObservation observationC;
    };
}