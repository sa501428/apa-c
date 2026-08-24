#include "apa_hic.h"

#include "straw.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using PairKey = uint64_t;
constexpr int32_t kCompression = 1000;

struct LoopRecord {
    std::string chrom1;
    std::string chrom2;
    int32_t chr1_index = 0;
    int32_t chr2_index = 0;
    int64_t start1 = 0;
    int64_t end1 = 0;
    int64_t start2 = 0;
    int64_t end2 = 0;
    int32_t center_bin1 = 0;
    int32_t center_bin2 = 0;
    int32_t window_start1 = 0;
    int32_t window_start2 = 0;
};

struct OutputPlan {
    OutputSpec spec;
    std::unordered_map<PairKey, std::vector<LoopRecord>> loops_by_pair;
};

struct OutputBuffers {
    std::vector<float> matrix;
    std::vector<double> row_sums;
    std::vector<double> col_sums;
};

struct PairWork {
    PairKey key = 0;
    chromosome chr1;
    chromosome chr2;
    int32_t num_bins1 = 0;
    int32_t num_bins2 = 0;
    std::vector<std::vector<LoopRecord>> loops_per_output;
    size_t total_loops = 0;
};

struct Bounds2DInfo {
    int32_t origin_start_x = 0;
    int32_t origin_start_y = 0;
    int32_t valid_start_x = 0;
    int32_t valid_end_x = 0;
    int32_t valid_start_y = 0;
    int32_t valid_end_y = 0;

    bool contains(int32_t bin_x, int32_t bin_y) const {
        return bin_x >= valid_start_x && bin_x < valid_end_x &&
               bin_y >= valid_start_y && bin_y < valid_end_y;
    }
};

class CompressedGridMapForBounds2D {
public:
    void build(const std::vector<LoopRecord>& loops,
               int32_t num_bins1,
               int32_t num_bins2,
               int matrix_width,
               std::vector<uint8_t>& all_row_indices,
               std::vector<uint8_t>& all_col_indices) {
        if (loops.empty()) {
            return;
        }

        grid_.reserve(loops.size() * 4);

        for (const auto& loop : loops) {
            const int32_t origin_start_x = loop.window_start1;
            const int32_t origin_start_y = loop.window_start2;
            const int32_t origin_end_x = origin_start_x + matrix_width;
            const int32_t origin_end_y = origin_start_y + matrix_width;

            const int32_t valid_start_x = std::max(0, origin_start_x);
            const int32_t valid_start_y = std::max(0, origin_start_y);
            const int32_t valid_end_x = std::min(num_bins1 + 1, origin_end_x);
            const int32_t valid_end_y = std::min(num_bins2 + 1, origin_end_y);

            if (valid_start_x >= valid_end_x || valid_start_y >= valid_end_y) {
                continue;
            }

            const Bounds2DInfo bounds{
                origin_start_x,
                origin_start_y,
                valid_start_x,
                valid_end_x,
                valid_start_y,
                valid_end_y,
            };

            const int32_t compressed_start_x = valid_start_x / kCompression;
            const int32_t compressed_start_y = valid_start_y / kCompression;
            const int32_t compressed_end_x = (valid_end_x - 1) / kCompression;
            const int32_t compressed_end_y = (valid_end_y - 1) / kCompression;

            for (int32_t row = compressed_start_x; row <= compressed_end_x; ++row) {
                for (int32_t col = compressed_start_y; col <= compressed_end_y; ++col) {
                    grid_[makeCompressedKey(row, col)].push_back(bounds);
                }
            }

            std::fill(all_row_indices.begin() + valid_start_x, all_row_indices.begin() + valid_end_x, 1U);
            std::fill(all_col_indices.begin() + valid_start_y, all_col_indices.begin() + valid_end_y, 1U);
        }
    }

    bool empty() const {
        return grid_.empty();
    }

    template <typename Fn>
    void forEachContainingBound(int32_t bin_x, int32_t bin_y, Fn&& fn) const {
        const auto it = grid_.find(makeCompressedKey(bin_x / kCompression, bin_y / kCompression));
        if (it == grid_.end()) {
            return;
        }

        for (const auto& bounds : it->second) {
            if (bounds.contains(bin_x, bin_y)) {
                fn(bounds);
            }
        }
    }

private:
    static PairKey makeCompressedKey(int32_t row, int32_t col) {
        return (static_cast<PairKey>(static_cast<uint32_t>(row)) << 32U) |
               static_cast<uint32_t>(col);
    }

    std::unordered_map<PairKey, std::vector<Bounds2DInfo>> grid_;
};

struct RegionsOfInterest {
    std::vector<CompressedGridMapForBounds2D> loop_maps;
    std::vector<uint8_t> all_row_indices;
    std::vector<uint8_t> all_col_indices;

    bool probablyContainsRecord(int32_t bin_x, int32_t bin_y) const {
        return bin_x >= 0 && bin_y >= 0 &&
               static_cast<size_t>(bin_x) < all_row_indices.size() &&
               static_cast<size_t>(bin_y) < all_col_indices.size() &&
               all_row_indices[static_cast<size_t>(bin_x)] != 0 &&
               all_col_indices[static_cast<size_t>(bin_y)] != 0;
    }
};

struct DistanceBounds {
    int32_t min_dist = 0;
    int32_t max_dist = std::numeric_limits<int32_t>::max();
};

struct ThreadLocalState {
    std::vector<OutputBuffers> outputs;
};

PairKey makePairKey(int32_t chr1_index, int32_t chr2_index) {
    return (static_cast<PairKey>(static_cast<uint32_t>(chr1_index)) << 32U) |
           static_cast<uint32_t>(chr2_index);
}

std::string toLower(std::string value) {
    for (char& ch : value) {
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    }
    return value;
}

std::string toUpper(std::string value) {
    for (char& ch : value) {
        ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
    }
    return value;
}

bool isAggNormLabel(const std::string& value) {
    const std::string lower = toLower(value);
    return lower.find("ag") != std::string::npos &&
           lower.find("norm") != std::string::npos;
}

bool endsWith(const std::string& value, const std::string& suffix) {
    return suffix.size() <= value.size() &&
           value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
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

std::unordered_map<std::string, chromosome> buildChromosomeMap(const std::vector<chromosome>& chromosomes) {
    std::unordered_map<std::string, chromosome> by_name;
    by_name.reserve(chromosomes.size());
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

void canonicalizeLoop(LoopRecord& loop, const std::unordered_map<std::string, chromosome>& chromosome_map) {
    const auto chr1_it = chromosome_map.find(loop.chrom1);
    const auto chr2_it = chromosome_map.find(loop.chrom2);
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
        std::swap(loop.start1, loop.start2);
        std::swap(loop.end1, loop.end2);
    }

    loop.chr1_index = chromosome_map.at(loop.chrom1).index;
    loop.chr2_index = chromosome_map.at(loop.chrom2).index;
}

bool loopPassesFilters(const LoopRecord& loop, const Options& options) {
    if (loop.chrom1 != loop.chrom2) {
        return options.include_inter;
    }

    const int32_t dist_bins = std::abs(loop.center_bin2 - loop.center_bin1);
    return dist_bins >= options.min_dist_bins && dist_bins <= options.max_dist_bins;
}

OutputPlan loadLoopsForOutput(const OutputSpec& spec,
                              const Options& options,
                              const std::unordered_map<std::string, chromosome>& chromosome_map) {
    std::ifstream in(spec.loop_file);
    if (!in) {
        throw std::runtime_error("Cannot open loop file: " + spec.loop_file);
    }

    OutputPlan plan;
    plan.spec = spec;

    std::string line;
    size_t line_number = 0;
    size_t kept_loops = 0;

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
        loop.window_start1 = loop.center_bin1 - options.window;
        loop.window_start2 = loop.center_bin2 - options.window;

        if (!loopPassesFilters(loop, options)) {
            continue;
        }

        plan.loops_by_pair[makePairKey(loop.chr1_index, loop.chr2_index)].push_back(loop);
        ++kept_loops;
    }

    if (kept_loops == 0) {
        throw std::runtime_error("No loops remained after filtering for " + spec.loop_file);
    }

    return plan;
}

std::vector<PairWork> buildPairWork(const std::vector<OutputPlan>& output_plans,
                                    const std::vector<chromosome>& chromosomes,
                                    int32_t resolution) {
    std::unordered_map<PairKey, size_t> pair_to_index;
    std::vector<PairWork> pair_work;

    for (size_t output_idx = 0; output_idx < output_plans.size(); ++output_idx) {
        for (const auto& entry : output_plans[output_idx].loops_by_pair) {
            const PairKey key = entry.first;
            const auto& loops = entry.second;
            if (loops.empty()) {
                continue;
            }

            auto it = pair_to_index.find(key);
            if (it == pair_to_index.end()) {
                const LoopRecord& first_loop = loops.front();
                PairWork work;
                work.key = key;
                work.chr1 = chromosomes[static_cast<size_t>(first_loop.chr1_index)];
                work.chr2 = chromosomes[static_cast<size_t>(first_loop.chr2_index)];
                work.num_bins1 = static_cast<int32_t>(work.chr1.length / resolution);
                work.num_bins2 = static_cast<int32_t>(work.chr2.length / resolution);
                work.loops_per_output.resize(output_plans.size());
                pair_work.push_back(std::move(work));
                it = pair_to_index.emplace(key, pair_work.size() - 1).first;
            }

            PairWork& work = pair_work[it->second];
            work.total_loops += loops.size();
            work.loops_per_output[output_idx] = loops;
        }
    }

    std::sort(pair_work.begin(), pair_work.end(), [](const PairWork& lhs, const PairWork& rhs) {
        if (lhs.total_loops != rhs.total_loops) {
            return lhs.total_loops > rhs.total_loops;
        }
        if (lhs.chr1.index != rhs.chr1.index) {
            return lhs.chr1.index < rhs.chr1.index;
        }
        return lhs.chr2.index < rhs.chr2.index;
    });

    return pair_work;
}

DistanceBounds calculateDistanceBounds(const PairWork& pair_work, int window) {
    DistanceBounds bounds;
    if (pair_work.chr1.index != pair_work.chr2.index) {
        return bounds;
    }

    int32_t min_dist = std::numeric_limits<int32_t>::max();
    int32_t max_dist = 0;

    for (const auto& loops : pair_work.loops_per_output) {
        for (const auto& loop : loops) {
            const int32_t dist = std::abs(loop.center_bin2 - loop.center_bin1);
            min_dist = std::min(min_dist, dist);
            max_dist = std::max(max_dist, dist);
        }
    }

    if (min_dist == std::numeric_limits<int32_t>::max()) {
        return bounds;
    }

    bounds.min_dist = std::max(0, min_dist - (2 * window));
    bounds.max_dist = max_dist + (2 * window);
    return bounds;
}

RegionsOfInterest buildRegionsOfInterest(const PairWork& pair_work, int matrix_width) {
    RegionsOfInterest roi;
    roi.all_row_indices.assign(static_cast<size_t>(pair_work.num_bins1) + 1, 0);
    roi.all_col_indices.assign(static_cast<size_t>(pair_work.num_bins2) + 1, 0);
    roi.loop_maps.resize(pair_work.loops_per_output.size());

    for (size_t output_idx = 0; output_idx < pair_work.loops_per_output.size(); ++output_idx) {
        roi.loop_maps[output_idx].build(pair_work.loops_per_output[output_idx],
                                        pair_work.num_bins1,
                                        pair_work.num_bins2,
                                        matrix_width,
                                        roi.all_row_indices,
                                        roi.all_col_indices);
    }

    return roi;
}

ThreadLocalState createThreadLocalState(size_t num_outputs, int matrix_width) {
    ThreadLocalState state;
    state.outputs.resize(num_outputs);

    const size_t matrix_size = static_cast<size_t>(matrix_width) * matrix_width;
    for (auto& output : state.outputs) {
        output.matrix.assign(matrix_size, 0.0f);
        output.row_sums.assign(matrix_width, 0.0);
        output.col_sums.assign(matrix_width, 0.0);
    }

    return state;
}

void addLocalNormSums(std::vector<double>& sums, const std::vector<double>& norm_vector, int32_t bin_start) {
    for (size_t i = 0; i < sums.size(); ++i) {
        const int64_t bin = static_cast<int64_t>(bin_start) + static_cast<int64_t>(i);
        if (bin < 0 || static_cast<size_t>(bin) >= norm_vector.size()) {
            continue;
        }

        const double value = norm_vector[static_cast<size_t>(bin)];
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

void applyAggregateNormalization(OutputBuffers& output, int matrix_width) {
    normalizeByAverage(output.row_sums);
    normalizeByAverage(output.col_sums);

    for (int r = 0; r < matrix_width; ++r) {
        for (int c = 0; c < matrix_width; ++c) {
            const double norm_value = output.row_sums[static_cast<size_t>(r)] *
                                      output.col_sums[static_cast<size_t>(c)];
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

void mergeOutputBuffers(std::vector<OutputBuffers>& global_outputs,
                        const std::vector<OutputBuffers>& local_outputs) {
    for (size_t output_idx = 0; output_idx < global_outputs.size(); ++output_idx) {
        auto& global = global_outputs[output_idx];
        const auto& local = local_outputs[output_idx];

        for (size_t i = 0; i < global.matrix.size(); ++i) {
            global.matrix[i] += local.matrix[i];
        }
        for (size_t i = 0; i < global.row_sums.size(); ++i) {
            global.row_sums[i] += local.row_sums[i];
            global.col_sums[i] += local.col_sums[i];
        }
    }
}

unsigned int getThreadCount(const Options& options, size_t num_pairs) {
    unsigned int configured_threads = options.threads > 0
        ? static_cast<unsigned int>(options.threads)
        : std::thread::hardware_concurrency();
    if (configured_threads == 0) {
        configured_threads = 1;
    }
    return std::max(1U, std::min(configured_threads, static_cast<unsigned int>(std::max<size_t>(1, num_pairs))));
}

void processPairWork(const PairWork& pair_work,
                     const Options& options,
                     ThreadLocalState& local_state,
                     std::mutex& io_mutex) {
    const int matrix_width = options.window * 2 + 1;
    const bool is_intra = pair_work.chr1.index == pair_work.chr2.index;

    {
        std::lock_guard<std::mutex> lock(io_mutex);
        std::cout << "Processing " << pair_work.chr1.name << " / " << pair_work.chr2.name
                  << " across " << pair_work.total_loops << " loops" << std::endl;
    }

    const RegionsOfInterest roi = buildRegionsOfInterest(pair_work, matrix_width);
    const DistanceBounds distance_bounds = calculateDistanceBounds(pair_work, options.window);

    const bool use_hic_normalization = !options.hic_normalization.empty();
    std::vector<double> custom_vc_norm_vector1;
    std::vector<double> custom_vc_norm_vector2;
    if (!use_hic_normalization) {
        custom_vc_norm_vector1.assign(static_cast<size_t>(pair_work.num_bins1) + 1, 0.0);
        custom_vc_norm_vector2.assign(static_cast<size_t>(pair_work.num_bins2) + 1, 0.0);
    }

    const auto process_record = [&](int32_t bin_x, int32_t bin_y, float counts) {
        if (!(counts > 0.0f) || !std::isfinite(counts) || bin_x < 0 || bin_y < 0) {
            return;
        }
        if (bin_x > pair_work.num_bins1 || bin_y > pair_work.num_bins2) {
            return;
        }

        if (!use_hic_normalization) {
            custom_vc_norm_vector1[static_cast<size_t>(bin_x)] += counts;
            custom_vc_norm_vector2[static_cast<size_t>(bin_y)] += counts;
        }

        if (is_intra) {
            const int32_t dist = std::abs(bin_y - bin_x);
            if (dist < distance_bounds.min_dist || dist > distance_bounds.max_dist) {
                return;
            }
        }
        if (!roi.probablyContainsRecord(bin_x, bin_y)) {
            return;
        }

        for (size_t output_idx = 0; output_idx < roi.loop_maps.size(); ++output_idx) {
            const auto& loop_map = roi.loop_maps[output_idx];
            if (loop_map.empty()) {
                continue;
            }

            auto& output_matrix = local_state.outputs[output_idx].matrix;
            loop_map.forEachContainingBound(bin_x, bin_y, [&](const Bounds2DInfo& bounds) {
                const int32_t relative_x = bin_x - bounds.origin_start_x;
                const int32_t relative_y = bin_y - bounds.origin_start_y;
                output_matrix[static_cast<size_t>(relative_x) * matrix_width + relative_y] += counts;
            });
        }
    };

    if (use_hic_normalization) {
        // strawStream returns genomic positions, whereas the raw block API returns bin indices.
        size_t normalized_record_count = 0;
        strawStream("observed",
                    options.hic_normalization,
                    options.hic_file,
                    pair_work.chr1.name,
                    pair_work.chr2.name,
                    "BP",
                    options.resolution,
                    [&](const contactRecord& record) {
            if (record.counts > 0.0f && std::isfinite(record.counts)) {
                ++normalized_record_count;
            }
            process_record(record.binX / options.resolution,
                           record.binY / options.resolution,
                           record.counts);
        });

        // straw reports a missing normalization by returning no records. Distinguish that
        // from an actually empty matrix so a misspelled or unavailable norm cannot silently
        // produce an all-zero APA result. This extra raw scan only runs on the error/empty path.
        if (normalized_record_count == 0) {
            bool raw_matrix_has_records = false;
            forEachRawObservedBlock(options.hic_file,
                                    pair_work.chr1.name,
                                    pair_work.chr2.name,
                                    options.resolution,
                                    [&](const std::vector<contactRecord>& block_records) {
                if (raw_matrix_has_records) {
                    return;
                }
                raw_matrix_has_records = std::any_of(
                    block_records.begin(), block_records.end(), [](const contactRecord& record) {
                        return record.counts > 0.0f && std::isfinite(record.counts);
                    });
            });
            if (raw_matrix_has_records) {
                throw std::runtime_error(
                    ".hic normalization " + options.hic_normalization +
                    " is unavailable or has no usable values for " + pair_work.chr1.name +
                    " / " + pair_work.chr2.name + " at resolution " +
                    std::to_string(options.resolution));
            }
        }
    } else {
        forEachRawObservedBlock(options.hic_file,
                                pair_work.chr1.name,
                                pair_work.chr2.name,
                                options.resolution,
                                [&](const std::vector<contactRecord>& block_records) {
            for (const auto& record : block_records) {
                process_record(record.binX, record.binY, record.counts);
            }
        });

        for (size_t output_idx = 0; output_idx < pair_work.loops_per_output.size(); ++output_idx) {
            auto& output = local_state.outputs[output_idx];
            for (const auto& loop : pair_work.loops_per_output[output_idx]) {
                addLocalNormSums(output.row_sums, custom_vc_norm_vector1, loop.window_start1);
                addLocalNormSums(output.col_sums, custom_vc_norm_vector2, loop.window_start2);
            }
        }
    }

    {
        std::lock_guard<std::mutex> lock(io_mutex);
        std::cout << "Completed " << pair_work.chr1.name << "_" << pair_work.chr2.name << std::endl;
    }
}

}  // namespace

void printUsage(std::ostream& out) {
    out << "Usage: apa_hic [--window W] [--min-dist BINS] [--max-dist BINS] [--threads N]\n"
        << "               [--include-inter] [-r RESOLUTION] [-k HIC_NORM]\n"
        << "               <input.hic> <loops1.bedpe> <output1.{npy,txt}> [<loops2.bedpe> <output2...> ...]\n"
        << "\n"
        << "Notes:\n"
        << "  --min-dist/--max-dist are interpreted in bins, matching APA3-style behavior.\n"
        << "  The default is locally calculated aggregate normalization.\n"
        << "  -k/--norm uses a normalization stored in the .hic file (for example KR or SCALE).\n"
        << "  --ag-norm explicitly selects the default aggregate normalization.\n"
        << "  Loop files are BEDPE-like; only the first 6 columns are required.\n";
}

Options parseOptions(int argc, char* argv[]) {
    Options options;
    std::vector<std::string> positional;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--ag-norm") {
            options.hic_normalization.clear();
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
        } else if (arg == "--threads") {
            options.threads = parseIntOption(nextValue(argc, argv, i, arg), arg);
        } else if (arg == "-k" || arg == "--norm") {
            const std::string norm_arg = nextValue(argc, argv, i, arg);
            if (norm_arg.empty()) {
                throw std::runtime_error("Normalization name cannot be empty");
            }
            options.hic_normalization = isAggNormLabel(norm_arg) ? "" : toUpper(norm_arg);
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
    if (options.threads < 0) {
        throw std::runtime_error("threads must be non-negative");
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
    const auto chromosome_map = buildChromosomeMap(chromosomes);
    const int matrix_width = options.window * 2 + 1;

    std::vector<OutputPlan> output_plans;
    output_plans.reserve(options.outputs.size());
    for (const auto& spec : options.outputs) {
        output_plans.push_back(loadLoopsForOutput(spec, options, chromosome_map));
    }

    std::vector<PairWork> pair_work = buildPairWork(output_plans, chromosomes, options.resolution);
    if (pair_work.empty()) {
        throw std::runtime_error("No chromosome pairs remained after loading loop files");
    }

    std::vector<OutputBuffers> global_outputs = createThreadLocalState(options.outputs.size(), matrix_width).outputs;

    std::atomic<size_t> next_pair_index{0};
    std::mutex io_mutex;
    std::mutex merge_mutex;
    std::mutex error_mutex;
    std::exception_ptr first_error;
    std::atomic<bool> stop_requested{false};

    const unsigned int num_threads = getThreadCount(options, pair_work.size());
    std::vector<std::thread> workers;
    workers.reserve(num_threads);

    for (unsigned int thread_index = 0; thread_index < num_threads; ++thread_index) {
        workers.emplace_back([&]() {
            ThreadLocalState local_state = createThreadLocalState(options.outputs.size(), matrix_width);
            try {
                while (!stop_requested.load(std::memory_order_relaxed)) {
                    const size_t pair_index = next_pair_index.fetch_add(1, std::memory_order_relaxed);
                    if (pair_index >= pair_work.size()) {
                        break;
                    }

                    processPairWork(pair_work[pair_index], options, local_state, io_mutex);
                }
            } catch (...) {
                stop_requested.store(true, std::memory_order_relaxed);
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error) {
                    first_error = std::current_exception();
                }
            }

            std::lock_guard<std::mutex> lock(merge_mutex);
            mergeOutputBuffers(global_outputs, local_state.outputs);
        });
    }

    for (auto& worker : workers) {
        worker.join();
    }

    if (first_error) {
        std::rethrow_exception(first_error);
    }

    for (size_t output_idx = 0; output_idx < global_outputs.size(); ++output_idx) {
        if (options.hic_normalization.empty()) {
            applyAggregateNormalization(global_outputs[output_idx], matrix_width);
        }
        std::cout << "Writing " << options.outputs[output_idx].output_file << std::endl;
        saveMatrix(options.outputs[output_idx].output_file, global_outputs[output_idx].matrix, matrix_width);
    }

    return 0;
}
