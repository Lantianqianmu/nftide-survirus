#include <cstring>
#include "ssw_cpp.h"

// Same scoring and aligned-query-span test as SurVirus filter.cpp.
extern "C" double aligned_fraction(const char* reference, const char* query) {
    if (!*query || !*reference) return 0;
    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;
    aligner.Align(query, reference, std::strlen(reference), filter, &alignment, 0);
    if (alignment.sw_score <= 0) return 0;
    // SSW reports inclusive endpoints. A 30-base exact match spans 30 bases.
    return double(alignment.query_end - alignment.query_begin + 1) / std::strlen(query);
}
