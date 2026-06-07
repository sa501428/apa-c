#include "apa_hic.h"

#include "../straw-master/C++/straw.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

using PairKey = std::pair<std::string, std::string>;

struct LoopRecord {
    std::string chrom1;
    std::string chrom2;
    int64_t start1 = 0;
    int64_t end1 = 0;
    int64_t start2 = 0;
    int64_t end2 = 0;
    int32_t center_bin1 = 0;
    int32_t center_bin2 = 0;
};

struct OutputState {
    OutputSpec spec;
    std::map<PairKey, std::vector<LoopRecord>> loops_by_pair;
    std::vector<float> matrix;
    std::vector<double> row_sums;
    std::vector<double> col_sums;
};

struct SparseMaskedContactMatrix {
    std::unordered_map<int32_t, std::unordered_map<int32_t, float>> data;

    void add(int32_t bin_x, int32_t bin_y, float counts) {
        data[bin_x][bin_y] = counts;
    }

    float get(int32_t bin_x, int32_t bin_y) const {
        auto row_it = data.find(bin_x);
        if (row_it == data.end()) {
            return 0.0f;
        }

        auto col_it = row_it->second.find(bin_y);
        if (col_it == row_it->second.end()) {
            return 0.0f;
        }

        return col_it->second;
    }

    void addLocalBoundedRegion(std::vector<float>& output, int matrix_width,
                               int32_t bin_x_start, int32_t bin_y_start) const {
        for (int r = 0; r < matrix_width; ++r) {
            for (int c = 0; c < matrix_width; ++c) {
                output[static_cast<size_t>(r) * matrix_width + c] +=
                    get(bin_x_start + r, bin_y_start + c);
            }
        }
    }
};

std::string toLower(std::string value) {
    for (char& ch : value) {
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    }
    return value;
}

bool isAggNormLabel(const std::string& value) {
    std::string lower = toLower(value);
    return lower.find("ag") != std::string::npos &&
           lower.find("norm") != std::string::npos;
}

bool endsWith(const std::string& value, const std::string& suffix) {
    if (suffix.size() > value.size()) {
        return false;
    }
    return value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

void requireValue(int argc, int index, const std::string& option) {
    if (index + 1 >= argc) {
        throw std::runtime_error("Missing value for option " + option);
    }
}

std::string nextValue(int argc, char* argv[], int& index, const std::string& option) {
    requireValue(argc, index, option);
    ++index;
    return argv[index];
}

int parseIntOption(const std::string& value, const std::string& option) {
    try {
        return std::stoi(value);
    } catch (const std::exception&) {
        throw std::runtime_error("Invalid integer for " + option + ": " + value);
    }
}

std::map<std::string, chromosome> buildChromosomeMap(const std::vector<chromosome>& chromosomes) {
    std::map<std::string, chromosome> by_name;
    for (const auto& chrom : chromosomes) {
        by_name[chrom.name] = chrom;
    }
    return by_name;
}

void validateResolutionAvailable(const Options& options) {
    const std::vector<int32_t> resolutions = getResolutionsForFile(options.hic_file);
    if (std::find(resolutions.begin(), resolutions.end(), options.resolution) == resolutions.end()) {
        std::ostringstream oss;
        oss << "Resolution " << options.resolution << " is not available in " << options.hic_file << ". Available:";
        for (int32_t resolution : resolutions) {
            oss << " " << resolution;
        }
        throw std::runtime_error(oss.str());
    }
}

LoopRecord parseLoopLine(const std::string& line, const std::string& file_name, size_t line_number) {
    std::istringstream iss(line);
    LoopRecord loop;
    if (!(iss >> loop.chrom1 >> loop.start1 >> loop.end1 >> loop.chrom2 >> loop.start2 >> loop.end2)) {
        throw std::runtime_error("Invalid BEDPE line in " + file_name + " at line " +
                                 std::to_string(line_number));
    }

    if (loop.end1 < loop.start1 || loop.end2 < loop.start2) {
        throw std::runtime_error("Invalid BEDPE interval in " + file_name + " at line " +
                                 std::to_string(line_number));
    }

    return loop;
}

void canonicalizeLoop(LoopRecord& loop, const std::map<std::string, chromosome>& chromosome_map) {
    auto chr1_it = chromosome_map.find(loop.chrom1);
    auto chr2_it = chromosome_map.find(loop.chrom2);
    if (chr1_it == chromosome_map.end() || chr2_it == chromosome_map.end()) {
        throw std::runtime_error("Loop references chromosome absent from .hic file: " +
                                 loop.chrom1 + " / " + loop.chrom2);
    }

    const int64_t mid1 = (loop.start1 + loop.end1) / 2;
    const int64_t mid2 = (loop.start2 + loop.end2) / 2;

    if (chr1_it->second.index > chr2_it->second.index) {
        std::swap(loop.chrom1, loop.chrom2);
        std::swap(loop.start1, loop.start2);
        std::swap(loop.end1, loop.end2);
    } else if (chr1_it->second.index == chr2_it->second.index &&
               (mid1 > mid2 || (mid1 == mid2 && loop.start1 > loop.start2))) {
        // For same-chromosome loops ensure the first anchor bins at or before the second,
        // matching the upper-triangle storage in the .hic file.
        std::swap(loop.start1, loop.start2);
        std::swap(loop.end1, loop.end2);
    }
}

bool loopPassesFilters(const LoopRecord& loop, const Options& options) {
    const bool is_inter = loop.chrom1 != loop.chrom2;
    if (is_inter) {
        return options.include_inter;
    }

    const int32_t dist_bins = std::abs(loop.center_bin2 - loop.center_bin1);
    return dist_bins >= options.min_dist_bins && dist_bins <= options.max_dist_bins;
}

std::vector<LoopRecord> loadLoopsForOutput(const OutputSpec& spec, const Options& options,
                                           const std::map<std::string, chromosome>& chromosome_map) {
    std::ifstream in(spec.loop_file);
    if (!in) {
        throw std::runtime_error("Cannot open loop file: " + spec.loop_file);
    }

    std::vector<LoopRecord> loops;
    std::string line;
    size_t line_number = 0;

    while (std::getline(in, line)) {
        ++line_number;
        if (line.empty() || line[0] == '#') {
            continue;
        }

        LoopRecord loop = parseLoopLine(line, spec.loop_file, line_number);
        if (chromosome_map.count(loop.chrom1) == 0 || chromosome_map.count(loop.chrom2) == 0) {
            continue;
        }

        canonicalizeLoop(loop, chromosome_map);

        loop.center_bin1 = static_cast<int32_t>(((loop.start1 + loop.end1) / 2) / options.resolution);
        loop.center_bin2 = static_cast<int32_t>(((loop.start2 + loop.end2) / 2) / options.resolution);

        if (loopPassesFilters(loop, options)) {
            loops.push_back(loop);
        }
    }

    if (loops.empty()) {
        throw std::runtime_error("No loops remained after filtering for " + spec.loop_file);
    }

    return loops;
}

void addLocalNormSums(std::vector<double>& sums, const std::vector<double>& norm_vector, int32_t bin_start) {
    for (size_t i = 0; i < sums.size(); ++i) {
        int64_t bin = static_cast<int64_t>(bin_start) + static_cast<int64_t>(i);
        if (bin < 0 || static_cast<size_t>(bin) >= norm_vector.size()) {
            continue;
        }

        double value = norm_vector[static_cast<size_t>(bin)];
        if (std::isfinite(value) && value > 0.0) {
            sums[i] += value;
        }
    }
}

double averageOf(const std::vector<double>& values) {
    if (values.empty()) {
        return 0.0;
    }

    double total = 0.0;
    for (double value : values) {
        total += value;
    }
    return total / static_cast<double>(values.size());
}

void normalizeByAverage(std::vector<double>& values) {
    const double average = averageOf(values);
    if (average <= 0.0) {
        return;
    }

    for (double& value : values) {
        value /= average;
    }
}

void applyAggregateNormalization(OutputState& output, int matrix_width) {
    normalizeByAverage(output.row_sums);
    normalizeByAverage(output.col_sums);

    for (int r = 0; r < matrix_width; ++r) {
        for (int c = 0; c < matrix_width; ++c) {
            const double norm_value = output.row_sums[r] * output.col_sums[c];
            float& cell = output.matrix[static_cast<size_t>(r) * matrix_width + c];
            if (norm_value > 0.0) {
                cell = static_cast<float>(cell / norm_value);
            } else {
                cell = 0.0f;
            }
        }
    }
}

void saveAsText(const std::string& file_name, const std::vector<float>& matrix, int matrix_width) {
    std::ofstream out(file_name);
    if (!out) {
        throw std::runtime_error("Cannot open output file: " + file_name);
    }

    out << std::fixed << std::setprecision(6);
    for (int r = 0; r < matrix_width; ++r) {
        for (int c = 0; c < matrix_width; ++c) {
            if (c > 0) {
                out << ' ';
            }
            out << matrix[static_cast<size_t>(r) * matrix_width + c];
        }
        out << '\n';
    }
}

void saveAsNpy(const std::string& file_name, const std::vector<float>& matrix, int matrix_width) {
    std::ofstream out(file_name, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Cannot open output file: " + file_name);
    }

    const std::string magic = "\x93NUMPY";
    const char major = 1;
    const char minor = 0;
    std::string header = "{'descr': '<f4', 'fortran_order': False, 'shape': (" +
                         std::to_string(matrix_width) + ", " + std::to_string(matrix_width) +
                         "), }";

    const size_t base_size = magic.size() + 2 + 2;
    const size_t remainder = (base_size + header.size() + 1) % 64;
    const size_t padding = (remainder == 0) ? 0 : 64 - remainder;
    header.append(padding, ' ');
    header.push_back('\n');

    const uint16_t header_len = static_cast<uint16_t>(header.size());

    out.write(magic.data(), static_cast<std::streamsize>(magic.size()));
    out.write(&major, 1);
    out.write(&minor, 1);
    out.write(reinterpret_cast<const char*>(&header_len), sizeof(header_len));
    out.write(header.data(), static_cast<std::streamsize>(header.size()));
    out.write(reinterpret_cast<const char*>(matrix.data()),
              static_cast<std::streamsize>(matrix.size() * sizeof(float)));
}

void saveMatrix(const std::string& file_name, const std::vector<float>& matrix, int matrix_width) {
    if (endsWith(file_name, ".npy")) {
        saveAsNpy(file_name, matrix, matrix_width);
    } else {
        saveAsText(file_name, matrix, matrix_width);
    }
}

SparseMaskedContactMatrix buildSparseMaskedMatrix(const std::vector<contactRecord>& records,
                                                  const std::unordered_set<int32_t>& relevant_bins,
                                                  int resolution, bool is_intra,
                                                  int min_dist_with_buffer, int max_dist_with_buffer,
                                                  std::vector<double>& coverage1,
                                                  std::vector<double>& coverage2) {
    SparseMaskedContactMatrix sparse;
    sparse.data.reserve(relevant_bins.size());

    for (const auto& record : records) {
        if (!(record.counts > 0.0f) || !std::isfinite(record.counts)) {
            continue;
        }

        const int32_t bin_x = record.binX / resolution;
        const int32_t bin_y = record.binY / resolution;

        if (bin_x >= 0) {
            if (static_cast<size_t>(bin_x) >= coverage1.size()) {
                coverage1.resize(static_cast<size_t>(bin_x) + 1, 0.0);
            }
            coverage1[static_cast<size_t>(bin_x)] += record.counts;
        }
        if (bin_y >= 0) {
            if (static_cast<size_t>(bin_y) >= coverage2.size()) {
                coverage2.resize(static_cast<size_t>(bin_y) + 1, 0.0);
            }
            coverage2[static_cast<size_t>(bin_y)] += record.counts;
        }

        if (is_intra) {
            const int32_t dist = std::abs(bin_y - bin_x);
            if (dist < min_dist_with_buffer || dist > max_dist_with_buffer) {
                continue;
            }
        }

        if (relevant_bins.count(bin_x) == 0 || relevant_bins.count(bin_y) == 0) {
            continue;
        }

        sparse.add(bin_x, bin_y, record.counts);
    }

    return sparse;
}

void populateOutputsForPair(const PairKey& pair_key, const Options& options,
                            std::vector<OutputState>& outputs) {
    const bool is_intra = pair_key.first == pair_key.second;
    const int matrix_width = options.window * 2 + 1;

    size_t total_loops = 0;
    int min_dist = std::numeric_limits<int>::max();
    int max_dist = 0;
    std::unordered_set<int32_t> relevant_bins;

    for (const auto& output : outputs) {
        auto it = output.loops_by_pair.find(pair_key);
        if (it == output.loops_by_pair.end()) {
            continue;
        }

        total_loops += it->second.size();
        for (const auto& loop : it->second) {
            const int32_t start1 = loop.center_bin1 - options.window;
            const int32_t start2 = loop.center_bin2 - options.window;
            for (int i = 0; i < matrix_width; ++i) {
                relevant_bins.insert(start1 + i);
                relevant_bins.insert(start2 + i);
            }

            if (is_intra) {
                const int32_t dist = std::abs(loop.center_bin2 - loop.center_bin1);
                min_dist = std::min(min_dist, dist);
                max_dist = std::max(max_dist, dist);
            }
        }
    }

    if (total_loops == 0) {
        return;
    }

    std::cout << "Processing " << pair_key.first << " / " << pair_key.second
              << " across " << total_loops << " loops" << std::endl;

    std::vector<contactRecord> records =
        straw("observed", "NONE", options.hic_file, pair_key.first, pair_key.second, "BP",
              options.resolution);

    int min_dist_with_buffer = 0;
    int max_dist_with_buffer = std::numeric_limits<int>::max();
    if (is_intra && min_dist != std::numeric_limits<int>::max()) {
        min_dist_with_buffer = std::max(0, min_dist - (4 * options.window));
        max_dist_with_buffer = max_dist + (4 * options.window);
    }

    std::vector<double> norm_vector1;
    std::vector<double> norm_vector2;

    const SparseMaskedContactMatrix sparse =
        buildSparseMaskedMatrix(records, relevant_bins, options.resolution, is_intra,
                                min_dist_with_buffer, max_dist_with_buffer,
                                norm_vector1, norm_vector2);

    for (auto& output : outputs) {
        auto it = output.loops_by_pair.find(pair_key);
        if (it == output.loops_by_pair.end()) {
            continue;
        }

        for (const auto& loop : it->second) {
            const int32_t start1 = loop.center_bin1 - options.window;
            const int32_t start2 = loop.center_bin2 - options.window;
            sparse.addLocalBoundedRegion(output.matrix, matrix_width, start1, start2);
            addLocalNormSums(output.row_sums, norm_vector1, start1);
            addLocalNormSums(output.col_sums, norm_vector2, start2);
        }
    }
}

}  // namespace

void printUsage(std::ostream& out) {
    out << "Usage: apa_hic [--window W] [--min-dist BINS] [--max-dist BINS]\n"
        << "               [--include-inter] [-r RESOLUTION]\n"
        << "               <input.hic> <loops1.bedpe> <output1.{npy,txt}> [<loops2.bedpe> <output2...> ...]\n"
        << "\n"
        << "Notes:\n"
        << "  --min-dist/--max-dist are interpreted in bins, matching APA3-style behavior.\n"
        << "  Aggregate normalization is always used.\n"
        << "  Coverage vectors are recomputed locally from raw observed counts.\n"
        << "  Loop files are BEDPE-like; only the first 6 columns are required.\n";
}

Options parseOptions(int argc, char* argv[]) {
    Options options;
    std::vector<std::string> positional;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--ag-norm") {
            continue;
        } else if (arg == "--include-inter") {
            options.include_inter = true;
        } else if (arg == "-r" || arg == "--res") {
            options.resolution = parseIntOption(nextValue(argc, argv, i, arg), arg);
        } else if (arg == "-w" || arg == "--window") {
            options.window = parseIntOption(nextValue(argc, argv, i, arg), arg);
        } else if (arg == "--min-dist") {
            options.min_dist_bins = parseIntOption(nextValue(argc, argv, i, arg), arg);
        } else if (arg == "--max-dist") {
            options.max_dist_bins = parseIntOption(nextValue(argc, argv, i, arg), arg);
        } else if (arg == "-k" || arg == "--norm") {
            const std::string norm_arg = nextValue(argc, argv, i, arg);
            if (!isAggNormLabel(norm_arg)) {
                throw std::runtime_error(
                    "Only aggregate normalization is supported; omit --norm or use --ag-norm.");
            }
        } else if (arg == "-h" || arg == "--help") {
            printUsage(std::cout);
            std::exit(0);
        } else if (!arg.empty() && arg[0] == '-') {
            throw std::runtime_error("Unknown option: " + arg);
        } else {
            positional.push_back(arg);
        }
    }

    if (options.resolution <= 0) {
        throw std::runtime_error("Resolution must be positive");
    }
    if (options.window <= 0) {
        throw std::runtime_error("Window must be positive");
    }

    if (options.min_dist_bins < 0) {
        options.min_dist_bins = 2 * options.window;
    }
    if (options.max_dist_bins < 0) {
        options.max_dist_bins = std::numeric_limits<int>::max();
    }
    if (options.max_dist_bins < options.min_dist_bins) {
        throw std::runtime_error("max-dist must be >= min-dist");
    }

    if (positional.size() < 3 || ((positional.size() - 1) % 2 != 0)) {
        std::ostringstream oss;
        printUsage(oss);
        throw std::runtime_error(oss.str());
    }

    options.hic_file = positional[0];
    for (size_t i = 1; i < positional.size(); i += 2) {
        options.outputs.push_back(OutputSpec{positional[i], positional[i + 1]});
    }

    return options;
}

int runAPA(const Options& options) {
    if (options.hic_file.empty()) {
        throw std::runtime_error("No .hic file provided");
    }

    validateResolutionAvailable(options);

    const std::vector<chromosome> chromosomes = getChromosomesForFile(options.hic_file);
    const std::map<std::string, chromosome> chromosome_map = buildChromosomeMap(chromosomes);
    const int matrix_width = options.window * 2 + 1;

    std::vector<OutputState> outputs;
    outputs.reserve(options.outputs.size());

    std::set<PairKey> all_pairs;
    for (const auto& spec : options.outputs) {
        OutputState state;
        state.spec = spec;
        state.matrix.assign(static_cast<size_t>(matrix_width) * matrix_width, 0.0f);
        state.row_sums.assign(matrix_width, 0.0);
        state.col_sums.assign(matrix_width, 0.0);

        const std::vector<LoopRecord> loops = loadLoopsForOutput(spec, options, chromosome_map);
        for (const auto& loop : loops) {
            const PairKey key(loop.chrom1, loop.chrom2);
            state.loops_by_pair[key].push_back(loop);
            all_pairs.insert(key);
        }

        outputs.push_back(std::move(state));
    }

    if (all_pairs.empty()) {
        throw std::runtime_error("No chromosome pairs remained after loading loop files");
    }

    for (const auto& pair_key : all_pairs) {
        populateOutputsForPair(pair_key, options, outputs);
    }

    for (auto& output : outputs) {
        applyAggregateNormalization(output, matrix_width);
        std::cout << "Writing " << output.spec.output_file << std::endl;
        saveMatrix(output.spec.output_file, output.matrix, matrix_width);
    }

    return 0;
}
