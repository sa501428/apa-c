#ifndef APA_HIC_H
#define APA_HIC_H

#include <iosfwd>
#include <string>
#include <vector>

struct OutputSpec {
    std::string loop_file;
    std::string output_file;
};

struct Options {
    bool include_inter = false;
    int resolution = 5000;
    int window = 10;
    int min_dist_bins = -1;
    int max_dist_bins = -1;
    std::string hic_file;
    std::vector<OutputSpec> outputs;
};

void printUsage(std::ostream& out);
Options parseOptions(int argc, char* argv[]);
int runAPA(const Options& options);

#endif
