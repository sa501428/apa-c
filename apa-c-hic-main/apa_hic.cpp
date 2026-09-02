#include "apa_hic.h"

#include "straw.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <charconv>
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
#include <zlib.h>

#ifdef HAVE_HTSLIB
#include <htslib/hts.h>
#include <htslib/sam.h>
#endif

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

    bool containsRow(int32_t bin_x) const {
        return bin_x >= 0 && static_cast<size_t>(bin_x) < all_row_indices.size() &&
               all_row_indices[static_cast<size_t>(bin_x)] != 0;
    }

    bool containsColumn(int32_t bin_y) const {
        return bin_y >= 0 && static_cast<size_t>(bin_y) < all_col_indices.size() &&
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

struct PairsLayout {
    int chrom1 = 1;
    int pos1 = 2;
    int chrom2 = 3;
    int pos2 = 4;
};

struct StreamPairState {
    const PairWork* work = nullptr;
    RegionsOfInterest roi;
    DistanceBounds distance_bounds;
    std::vector<double> vc1;
    std::vector<double> vc2;
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

std::string stripGzipSuffix(std::string value) {
    if (endsWith(toLower(value), ".gz")) {
        value.resize(value.size() - 3);
    }
    return value;
}

InputFormat detectInputFormat(const std::string& file_name) {
    const std::string lower = toLower(stripGzipSuffix(file_name));
    if (endsWith(lower, ".hic")) {
        return InputFormat::Hic;
    }
    if (endsWith(lower, ".pairs")) {
        return InputFormat::Pairs;
    }
    if (endsWith(lower, ".bam")) {
        return InputFormat::Bam;
    }
    throw std::runtime_error(
        "Cannot determine input format for " + file_name +
        "; use --input-format hic, pairs, or bam");
}

InputFormat parseInputFormat(const std::string& value) {
    const std::string lower = toLower(value);
    if (lower == "auto") return InputFormat::Auto;
    if (lower == "hic") return InputFormat::Hic;
    if (lower == "pairs") return InputFormat::Pairs;
    if (lower == "bam") return InputFormat::Bam;
    throw std::runtime_error("Unknown input format: " + value);
}

const char* inputFormatName(InputFormat format) {
    switch (format) {
        case InputFormat::Hic: return ".hic";
        case InputFormat::Pairs: return "pairs";
        case InputFormat::Bam: return "BAM";
        case InputFormat::Auto: break;
    }
    return "input";
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

class TextLineReader {
public:
    explicit TextLineReader(const std::string& file_name) : file_name_(file_name) {
        file_ = gzopen(file_name.c_str(), "rb");
        if (file_ == nullptr) {
            throw std::runtime_error("Cannot open pairs file: " + file_name);
        }
    }

    ~TextLineReader() {
        if (file_ != nullptr) {
            gzclose(file_);
        }
    }

    bool next(std::string& line) {
        line.clear();
        while (true) {
            char* result = gzgets(file_, buffer_.data(), static_cast<int>(buffer_.size()));
            if (result == nullptr) {
                if (gzeof(file_)) {
                    return !line.empty();
                }
                int error_number = Z_OK;
                const char* message = gzerror(file_, &error_number);
                throw std::runtime_error("Error reading " + file_name_ + ": " +
                                         (message == nullptr ? "unknown zlib error" : message));
            }

            line.append(buffer_.data());
            if (!line.empty() && line.back() == '\n') {
                line.pop_back();
                if (!line.empty() && line.back() == '\r') {
                    line.pop_back();
                }
                return true;
            }
        }
    }

private:
    std::string file_name_;
    gzFile file_ = nullptr;
    std::array<char, 64 * 1024> buffer_{};
};

std::vector<std::string_view> splitFields(const std::string& line) {
    std::vector<std::string_view> fields;
    const char* cursor = line.data();
    const char* end = cursor + line.size();
    while (cursor < end) {
        while (cursor < end && (*cursor == ' ' || *cursor == '\t')) {
            ++cursor;
        }
        if (cursor == end) {
            break;
        }
        const char* start = cursor;
        while (cursor < end && *cursor != ' ' && *cursor != '\t') {
            ++cursor;
        }
        fields.emplace_back(start, static_cast<size_t>(cursor - start));
    }
    return fields;
}

int64_t parseCoordinate(std::string_view field, const std::string& file_name, size_t line_number) {
    int64_t value = 0;
    const auto parsed = std::from_chars(field.data(), field.data() + field.size(), value);
    if (parsed.ec != std::errc() || parsed.ptr != field.data() + field.size()) {
        throw std::runtime_error("Invalid coordinate in " + file_name + " at line " +
                                 std::to_string(line_number));
    }
    return value;
}

bool updatePairsLayout(const std::string& line, PairsLayout& layout) {
    constexpr std::string_view prefix = "#columns:";
    if (line.size() < prefix.size() ||
        std::string_view(line.data(), prefix.size()) != prefix) {
        return false;
    }

    const std::string columns_text = line.substr(prefix.size());
    const auto columns = splitFields(columns_text);
    auto find_column = [&](std::initializer_list<std::string_view> names) {
        for (size_t i = 0; i < columns.size(); ++i) {
            for (std::string_view name : names) {
                if (columns[i] == name) {
                    return static_cast<int>(i);
                }
            }
        }
        return -1;
    };

    const PairsLayout declared{
        find_column({"chrom1", "chr1"}),
        find_column({"pos1"}),
        find_column({"chrom2", "chr2"}),
        find_column({"pos2"}),
    };
    if (declared.chrom1 < 0 || declared.pos1 < 0 ||
        declared.chrom2 < 0 || declared.pos2 < 0) {
        throw std::runtime_error(
            "pairs #columns header must contain chrom1, pos1, chrom2, and pos2");
    }
    layout = declared;
    return true;
}

bool parsePairsDataLine(const std::string& line,
                        const PairsLayout& layout,
                        const std::string& file_name,
                        size_t line_number,
                        std::string_view& chrom1,
                        int64_t& pos1,
                        std::string_view& chrom2,
                        int64_t& pos2) {
    if (line.empty() || line[0] == '#') {
        return false;
    }
    const int max_column = std::max({layout.chrom1, layout.pos1, layout.chrom2, layout.pos2});
    constexpr size_t max_supported_columns = 128;
    if (max_column < 0 || static_cast<size_t>(max_column) >= max_supported_columns) {
        throw std::runtime_error("Required pairs column is beyond the supported 128-column limit");
    }

    // Contact lines dominate runtime. Split only through the last required
    // column into stack-backed views, avoiding per-contact vector allocations.
    const std::array<int, 4> required_columns{
        layout.chrom1, layout.pos1, layout.chrom2, layout.pos2};
    std::array<std::string_view, 4> required_fields{};
    size_t field_count = 0;
    const char* cursor = line.data();
    const char* end = cursor + line.size();
    while (cursor < end && field_count <= static_cast<size_t>(max_column)) {
        while (cursor < end && (*cursor == ' ' || *cursor == '\t')) ++cursor;
        if (cursor == end) break;
        const char* start = cursor;
        while (cursor < end && *cursor != ' ' && *cursor != '\t') ++cursor;
        const std::string_view field(start, static_cast<size_t>(cursor - start));
        for (size_t required = 0; required < required_columns.size(); ++required) {
            if (required_columns[required] == static_cast<int>(field_count)) {
                required_fields[required] = field;
            }
        }
        ++field_count;
    }
    if (field_count <= static_cast<size_t>(max_column)) {
        throw std::runtime_error("Invalid pairs line in " + file_name + " at line " +
                                 std::to_string(line_number));
    }
    chrom1 = required_fields[0];
    pos1 = parseCoordinate(required_fields[1], file_name, line_number);
    chrom2 = required_fields[2];
    pos2 = parseCoordinate(required_fields[3], file_name, line_number);
    // 4DN/DCIC pairs coordinates are 1-based. Convert to the zero-based
    // coordinates used by BAM and by the binning logic below.
    pos1 = std::max<int64_t>(0, pos1 - 1);
    pos2 = std::max<int64_t>(0, pos2 - 1);
    return true;
}

std::vector<chromosome> getChromosomesForPairs(const std::string& file_name) {
    TextLineReader reader(file_name);
    PairsLayout layout;
    std::vector<chromosome> chromosomes;
    std::unordered_map<std::string, size_t> indices;

    auto add_chromosome = [&](std::string_view name, int64_t length) {
        if (name.empty() || name == "!") {
            return;
        }
        const std::string owned_name(name);
        const auto existing = indices.find(owned_name);
        if (existing == indices.end()) {
            const size_t index = chromosomes.size();
            indices.emplace(owned_name, index);
            chromosomes.push_back(chromosome{owned_name, static_cast<int32_t>(index), length});
        } else {
            chromosomes[existing->second].length =
                std::max(chromosomes[existing->second].length, length);
        }
    };

    bool has_declared_sizes = false;
    std::string line;
    size_t line_number = 0;
    while (reader.next(line)) {
        ++line_number;
        if (updatePairsLayout(line, layout)) {
            continue;
        }
        constexpr std::string_view chromsize_prefix = "#chromsize:";
        if (line.size() >= chromsize_prefix.size() &&
            std::string_view(line.data(), chromsize_prefix.size()) == chromsize_prefix) {
            const auto fields = splitFields(line.substr(chromsize_prefix.size()));
            if (fields.size() < 2) {
                throw std::runtime_error("Invalid #chromsize header in " + file_name +
                                         " at line " + std::to_string(line_number));
            }
            const int64_t length = parseCoordinate(fields[1], file_name, line_number);
            if (length <= 0) {
                throw std::runtime_error("Chromosome size must be positive in " + file_name +
                                         " at line " + std::to_string(line_number));
            }
            add_chromosome(fields[0], length);
            has_declared_sizes = true;
            continue;
        }
        if (line.empty() || line[0] == '#') {
            continue;
        }
        if (has_declared_sizes) {
            // A standard pairs header declares the chromosome dictionary. Avoid
            // an unnecessary full discovery scan of large inputs when it does.
            break;
        }

        std::string_view chrom1;
        std::string_view chrom2;
        int64_t pos1 = 0;
        int64_t pos2 = 0;
        if (parsePairsDataLine(line, layout, file_name, line_number,
                               chrom1, pos1, chrom2, pos2)) {
            add_chromosome(chrom1, pos1 + 1);
            add_chromosome(chrom2, pos2 + 1);
        }
    }

    if (chromosomes.empty()) {
        throw std::runtime_error("No chromosomes found in pairs file: " + file_name);
    }
    if (!has_declared_sizes) {
        std::cerr << "Warning: " << file_name
                  << " has no #chromsize header; chromosome lengths were inferred from contacts"
                  << std::endl;
    }
    return chromosomes;
}

template <typename Fn>
void forEachPairsContact(const std::string& file_name,
                         const std::unordered_map<std::string, chromosome>& chromosome_map,
                         Fn&& fn) {
    TextLineReader reader(file_name);
    PairsLayout layout;
    std::string cached_chrom1_name;
    std::string cached_chrom2_name;
    int32_t cached_chrom1_index = -1;
    int32_t cached_chrom2_index = -1;
    auto resolve_chromosome = [&](std::string_view name,
                                  std::string& cached_name,
                                  int32_t& cached_index) {
        if (name == cached_name) {
            return cached_index;
        }
        cached_name.assign(name.data(), name.size());
        const auto found = chromosome_map.find(cached_name);
        cached_index = found == chromosome_map.end() ? -1 : found->second.index;
        return cached_index;
    };
    std::string line;
    size_t line_number = 0;
    while (reader.next(line)) {
        ++line_number;
        if (updatePairsLayout(line, layout) || line.empty() || line[0] == '#') {
            continue;
        }
        std::string_view chrom1_name;
        std::string_view chrom2_name;
        int64_t pos1 = 0;
        int64_t pos2 = 0;
        if (!parsePairsDataLine(line, layout, file_name, line_number,
                                chrom1_name, pos1, chrom2_name, pos2)) {
            continue;
        }
        const int32_t chr1_index = resolve_chromosome(
            chrom1_name, cached_chrom1_name, cached_chrom1_index);
        const int32_t chr2_index = resolve_chromosome(
            chrom2_name, cached_chrom2_name, cached_chrom2_index);
        if (chr1_index < 0 || chr2_index < 0) {
            continue;
        }
        fn(chr1_index, pos1, chr2_index, pos2);
    }
}

std::vector<chromosome> getChromosomesForBam(const std::string& file_name) {
#ifdef HAVE_HTSLIB
    samFile* file = sam_open(file_name.c_str(), "rb");
    if (file == nullptr) {
        throw std::runtime_error("Cannot open BAM file: " + file_name);
    }
    sam_hdr_t* header = sam_hdr_read(file);
    if (header == nullptr) {
        sam_close(file);
        throw std::runtime_error("Cannot read BAM header: " + file_name);
    }
    std::vector<chromosome> chromosomes;
    chromosomes.reserve(static_cast<size_t>(header->n_targets));
    for (int32_t i = 0; i < header->n_targets; ++i) {
        chromosomes.push_back(chromosome{
            header->target_name[i], i, static_cast<int64_t>(header->target_len[i])});
    }
    sam_hdr_destroy(header);
    sam_close(file);
    if (chromosomes.empty()) {
        throw std::runtime_error("BAM header has no reference sequences: " + file_name);
    }
    return chromosomes;
#else
    (void)file_name;
    throw std::runtime_error(
        "BAM input is unavailable because apa_hic was built without htslib; install htslib and rebuild");
#endif
}

template <typename Fn>
void forEachBamContact(const Options& options, Fn&& fn) {
#ifdef HAVE_HTSLIB
    samFile* file = sam_open(options.input_file.c_str(), "rb");
    if (file == nullptr) {
        throw std::runtime_error("Cannot open BAM file: " + options.input_file);
    }
    if (options.threads > 1) {
        hts_set_threads(file, options.threads);
    }
    sam_hdr_t* header = sam_hdr_read(file);
    bam1_t* record = bam_init1();
    if (header == nullptr || record == nullptr) {
        if (record != nullptr) bam_destroy1(record);
        if (header != nullptr) sam_hdr_destroy(header);
        sam_close(file);
        throw std::runtime_error("Cannot initialize BAM reader for " + options.input_file);
    }

    int status = 0;
    while ((status = sam_read1(file, header, record)) >= 0) {
        const uint16_t flag = record->core.flag;
        const uint16_t excluded = BAM_FUNMAP | BAM_FMUNMAP | BAM_FSECONDARY |
                                  BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FDUP;
        if ((flag & excluded) != 0 || (flag & BAM_FPAIRED) == 0 ||
            (flag & BAM_FREAD1) == 0 || record->core.tid < 0 || record->core.mtid < 0 ||
            record->core.pos < 0 || record->core.mpos < 0) {
            continue;
        }
        fn(record->core.tid,
           static_cast<int64_t>(record->core.pos),
           record->core.mtid,
           static_cast<int64_t>(record->core.mpos));
    }
    bam_destroy1(record);
    sam_hdr_destroy(header);
    sam_close(file);
    if (status < -1) {
        throw std::runtime_error("Error reading BAM records from " + options.input_file);
    }
#else
    (void)options;
    (void)fn;
    throw std::runtime_error(
        "BAM input is unavailable because apa_hic was built without htslib; install htslib and rebuild");
#endif
}

void validateResolutionAvailable(const Options& options) {
    const std::vector<int32_t> resolutions = getResolutionsForFile(options.input_file);
    if (std::find(resolutions.begin(), resolutions.end(), options.resolution) == resolutions.end()) {
        std::ostringstream oss;
        oss << "Resolution " << options.resolution << " is not available in " << options.input_file << ". Available:";
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

void processHicPairWork(const PairWork& pair_work,
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

        const bool row_is_needed = roi.containsRow(bin_x);
        const bool column_is_needed = roi.containsColumn(bin_y);
        if (!use_hic_normalization && row_is_needed) {
            custom_vc_norm_vector1[static_cast<size_t>(bin_x)] += counts;
        }
        if (!use_hic_normalization && column_is_needed) {
            custom_vc_norm_vector2[static_cast<size_t>(bin_y)] += counts;
        }

        // A one-axis overlap is still needed for the corresponding VC value,
        // but a matrix contribution requires both axes to overlap an ROI.
        if (!row_is_needed || !column_is_needed) {
            return;
        }

        if (is_intra) {
            const int32_t dist = std::abs(bin_y - bin_x);
            if (dist < distance_bounds.min_dist || dist > distance_bounds.max_dist) {
                return;
            }
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
                    options.input_file,
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
            forEachRawObservedBlock(options.input_file,
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
        forEachRawObservedBlock(options.input_file,
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

void processStreamInput(const std::vector<PairWork>& pair_work,
                        const std::unordered_map<std::string, chromosome>& chromosome_map,
                        const Options& options,
                        std::vector<OutputBuffers>& outputs) {
    const int matrix_width = options.window * 2 + 1;
    std::vector<StreamPairState> states;
    states.reserve(pair_work.size());
    std::unordered_map<PairKey, size_t> state_by_pair;
    state_by_pair.reserve(pair_work.size());

    for (const auto& work : pair_work) {
        StreamPairState state;
        state.work = &work;
        state.roi = buildRegionsOfInterest(work, matrix_width);
        state.distance_bounds = calculateDistanceBounds(work, options.window);
        state.vc1.assign(static_cast<size_t>(work.num_bins1) + 1, 0.0);
        state.vc2.assign(static_cast<size_t>(work.num_bins2) + 1, 0.0);
        state_by_pair.emplace(work.key, states.size());
        states.push_back(std::move(state));
    }

    uint64_t input_contacts = 0;
    uint64_t relevant_contacts = 0;
    const auto process_contact = [&](int32_t chr1_index,
                                     int64_t pos1,
                                     int32_t chr2_index,
                                     int64_t pos2) {
        ++input_contacts;
        if (chr1_index > chr2_index ||
            (chr1_index == chr2_index && pos1 > pos2)) {
            std::swap(chr1_index, chr2_index);
            std::swap(pos1, pos2);
        }
        const auto state_it = state_by_pair.find(makePairKey(chr1_index, chr2_index));
        if (state_it == state_by_pair.end()) {
            return;
        }

        StreamPairState& state = states[state_it->second];
        const PairWork& work = *state.work;
        const int64_t bin_x64 = pos1 / options.resolution;
        const int64_t bin_y64 = pos2 / options.resolution;
        if (bin_x64 < 0 || bin_y64 < 0 ||
            bin_x64 > work.num_bins1 || bin_y64 > work.num_bins2) {
            return;
        }
        const int32_t bin_x = static_cast<int32_t>(bin_x64);
        const int32_t bin_y = static_cast<int32_t>(bin_y64);
        const bool row_is_needed = state.roi.containsRow(bin_x);
        const bool column_is_needed = state.roi.containsColumn(bin_y);
        if (!row_is_needed && !column_is_needed) {
            return;
        }
        ++relevant_contacts;
        if (row_is_needed) {
            state.vc1[static_cast<size_t>(bin_x)] += 1.0;
        }
        if (column_is_needed) {
            state.vc2[static_cast<size_t>(bin_y)] += 1.0;
        }

        // Contacts with only one relevant axis have now served their VC
        // purpose and can be discarded before distance and exact-ROI checks.
        if (!row_is_needed || !column_is_needed) {
            return;
        }

        if (chr1_index == chr2_index) {
            const int32_t distance = std::abs(bin_y - bin_x);
            if (distance < state.distance_bounds.min_dist ||
                distance > state.distance_bounds.max_dist) {
                return;
            }
        }
        for (size_t output_idx = 0; output_idx < state.roi.loop_maps.size(); ++output_idx) {
            const auto& loop_map = state.roi.loop_maps[output_idx];
            if (loop_map.empty()) {
                continue;
            }
            auto& matrix = outputs[output_idx].matrix;
            loop_map.forEachContainingBound(bin_x, bin_y, [&](const Bounds2DInfo& bounds) {
                const int32_t relative_x = bin_x - bounds.origin_start_x;
                const int32_t relative_y = bin_y - bounds.origin_start_y;
                matrix[static_cast<size_t>(relative_x) * matrix_width + relative_y] += 1.0f;
            });
        }
    };

    std::cout << "Scanning " << inputFormatName(options.input_format) << " contacts from "
              << options.input_file << std::endl;
    if (options.input_format == InputFormat::Pairs) {
        forEachPairsContact(options.input_file, chromosome_map, process_contact);
    } else if (options.input_format == InputFormat::Bam) {
        forEachBamContact(options, process_contact);
    } else {
        throw std::runtime_error("Internal error: stream processing requested for non-stream input");
    }

    for (auto& state : states) {
        const PairWork& work = *state.work;
        for (size_t output_idx = 0; output_idx < work.loops_per_output.size(); ++output_idx) {
            auto& output = outputs[output_idx];
            for (const auto& loop : work.loops_per_output[output_idx]) {
                addLocalNormSums(output.row_sums, state.vc1, loop.window_start1);
                addLocalNormSums(output.col_sums, state.vc2, loop.window_start2);
            }
        }
    }

    std::cout << "Read " << input_contacts << " usable contacts; " << relevant_contacts
              << " overlapped at least one retained ROI axis" << std::endl;
}

}  // namespace

void printUsage(std::ostream& out) {
    out << "Usage: apa_hic [--window W] [--min-dist BINS] [--max-dist BINS] [--threads N]\n"
        << "               [--include-inter] [-r RESOLUTION] [-k HIC_NORM]\n"
        << "               [--input-format auto|hic|pairs|bam]\n"
        << "               <input.{hic,pairs,pairs.gz,bam}> <loops1.bedpe> <output1.{npy,txt}>\n"
        << "               [<loops2.bedpe> <output2...> ...]\n"
        << "\n"
        << "Notes:\n"
        << "  --min-dist/--max-dist are interpreted in bins, matching APA3-style behavior.\n"
        << "  The default is locally calculated aggregate normalization.\n"
        << "  pairs and BAM inputs support aggregate normalization only.\n"
        << "  pairs coordinates follow the 1-based 4DN/DCIC convention.\n"
        << "  BAM input uses primary, non-duplicate read1 records and their mate positions.\n"
        << "  -k/--norm uses a normalization stored in a .hic file (for example KR or SCALE).\n"
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
        } else if (arg == "--input-format") {
            options.input_format = parseInputFormat(nextValue(argc, argv, i, arg));
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

    options.input_file = positional[0];
    if (options.input_format == InputFormat::Auto) {
        options.input_format = detectInputFormat(options.input_file);
    }
    for (size_t i = 1; i < positional.size(); i += 2) {
        options.outputs.push_back(OutputSpec{positional[i], positional[i + 1]});
    }

    return options;
}

int runAPA(const Options& options) {
    if (options.input_file.empty()) {
        throw std::runtime_error("No contact input file provided");
    }

    if (options.input_format != InputFormat::Hic && !options.hic_normalization.empty()) {
        throw std::runtime_error(
            "Stored normalization " + options.hic_normalization +
            " can only be used with .hic input; pairs and BAM support aggregate normalization only");
    }

    std::vector<chromosome> chromosomes;
    if (options.input_format == InputFormat::Hic) {
        validateResolutionAvailable(options);
        chromosomes = getChromosomesForFile(options.input_file);
    } else if (options.input_format == InputFormat::Pairs) {
        chromosomes = getChromosomesForPairs(options.input_file);
    } else if (options.input_format == InputFormat::Bam) {
        chromosomes = getChromosomesForBam(options.input_file);
    } else {
        throw std::runtime_error("Input format was not resolved");
    }
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

    if (options.input_format == InputFormat::Hic) {
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

                        processHicPairWork(pair_work[pair_index], options, local_state, io_mutex);
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
    } else {
        processStreamInput(pair_work, chromosome_map, options, global_outputs);
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
