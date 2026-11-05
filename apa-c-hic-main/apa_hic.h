#ifndef APA_HIC_H
#define APA_HIC_H

#include <iosfwd>
#include <string>
#include <vector>

struct OutputSpec {
    std::string loop_file;
    std::string output_file;
};

enum class InputFormat {
    Auto,
    Hic,
    Pairs,
    Bam,
};

struct Options {
    bool include_inter = false;
    int resolution = 5000;
    int window = 10;
    int min_dist_bins = -1;
    int max_dist_bins = -1;
    int threads = 0;
    // Empty means the default, locally calculated aggregate normalization.
    // Otherwise this is a normalization stored in the .hic file (for example KR or SCALE).
    std::string hic_normalization;
    InputFormat input_format = InputFormat::Auto;
    std::string input_file;
    std::vector<OutputSpec> outputs;
};

void printUsage(std::ostream& out);
Options parseOptions(int argc, char* argv[]);
int runAPA(const Options& options);

#endif
