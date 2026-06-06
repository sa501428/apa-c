#ifndef HIC_SLICE_H
#define HIC_SLICE_H

#include <string>
#include <map>
#include <vector>
#include "straw.h"

// Magic string to identify file format
const std::string HICSLICE_MAGIC = "HICSLICE";

// Header structure for the file
struct HicSliceHeader {
    int32_t resolution;
    int32_t numChromosomes;
    std::map<std::string, int16_t> chromosomeKeys;
};

// Record structure for storage
struct CompressedContactRecord {
    int16_t chr1Key;
    int32_t binX;
    int16_t chr2Key;
    int32_t binY;
    float value;
};

// Filter type for contact records
enum class ContactFilter {
    ALL,           // Default: save all records
    INTRA_SHORT,   // Only intra-chromosomal records within 5Mb
    INTRA_LONG,    // Only intra-chromosomal records beyond 5Mb
    INTER,         // Only inter-chromosomal records
    INTRA          // All intra-chromosomal records
};

void dumpGenomeWideDataAtResolution(const std::string& matrixType, 
                                  const std::string& norm, 
                                  const std::string& filePath, 
                                  const std::string& unit, 
                                  int32_t resolution, 
                                  const std::string& outputPath,
                                  bool compressed = true,
                                  ContactFilter filter = ContactFilter::ALL);

#endif 