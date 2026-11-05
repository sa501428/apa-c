// Standalone ATA implementation for BigWig + BED input.
//
// Build example:
//   ./ata-c-bigwig/build.sh (from the repository root)
//
// Usage:
//   ./ata_bigwig [--window N] [--res N] [--mode java|coverage]
//                [--normalization all|overlap] signal.bw peaks.bed out_prefix
//
// The optional "java" mode mirrors src/cli/clt/anchor/ATA.java:
//   - center each BED interval at (start + end) / 2
//   - query midpoint +/- window * resolution
//   - place each positive BigWig interval value into the output bin determined
//     by the interval start
//   - normalize by the number of BED intervals with at least one BigWig hit
//
// The default "coverage" mode is often more natural for raw BigWig data:
//   - distribute value * covered_bases into every overlapped output bin
//   - divide by bin width and then by the number of contributing BED intervals
//
// BED6 strand is respected in both modes. Minus-strand profiles are reversed
// around the interval midpoint before aggregation so all features are oriented
// in the plus-strand direction. Missing, '.', and '+' strands are unchanged.

#include <bigWig.h>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

struct Region {
    std::string chrom;
    int64_t start = 0;
    int64_t end = 0;
    char strand = '.';
};

struct Args {
    std::string bigwig;
    std::string bed;
    std::string out_prefix;
    int64_t window = 100;
    int64_t resolution = 1;
    bool coverage_mode = true;
    bool normalize_all = true;
};

static void usage(const char *argv0) {
    std::cerr
        << "Usage: " << argv0
        << " [--window N] [--res N] [--mode java|coverage]"
        << " [--normalization all|overlap] signal.bw peaks.bed out_prefix\n"
        << "Defaults: --window 100 --res 1 --mode coverage --normalization all\n";
}

static int64_t parse_positive_i64(const std::string &value, std::string_view name) {
    size_t idx = 0;
    long long parsed = std::stoll(value, &idx);
    if (idx != value.size() || parsed < 1) {
        throw std::runtime_error(std::string(name) + " must be a positive integer");
    }
    return parsed;
}

static Args parse_args(int argc, char **argv) {
    Args args;
    std::vector<std::string> positional;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--window") {
            if (++i >= argc) throw std::runtime_error("--window requires a value");
            args.window = parse_positive_i64(argv[i], "--window");
        } else if (arg == "--res") {
            if (++i >= argc) throw std::runtime_error("--res requires a value");
            args.resolution = parse_positive_i64(argv[i], "--res");
        } else if (arg == "--mode") {
            if (++i >= argc) throw std::runtime_error("--mode requires java or coverage");
            std::string mode = argv[i];
            if (mode == "java") {
                args.coverage_mode = false;
            } else if (mode == "coverage") {
                args.coverage_mode = true;
            } else {
                throw std::runtime_error("--mode requires java or coverage");
            }
        } else if (arg == "--normalization") {
            if (++i >= argc) throw std::runtime_error("--normalization requires all or overlap");
            std::string normalization = argv[i];
            if (normalization == "all") {
                args.normalize_all = true;
            } else if (normalization == "overlap") {
                args.normalize_all = false;
            } else {
                throw std::runtime_error("--normalization requires all or overlap");
            }
        } else if (arg == "-h" || arg == "--help") {
            usage(argv[0]);
            std::exit(0);
        } else {
            positional.push_back(std::move(arg));
        }
    }

    if (positional.size() != 3) {
        usage(argv[0]);
        throw std::runtime_error("expected signal.bw peaks.bed out_prefix");
    }

    args.bigwig = positional[0];
    args.bed = positional[1];
    args.out_prefix = positional[2];
    return args;
}

static bool starts_with(std::string_view s, std::string_view prefix) {
    return s.size() >= prefix.size() && s.substr(0, prefix.size()) == prefix;
}

static std::vector<std::string_view> split_ws(std::string_view line) {
    std::vector<std::string_view> fields;
    size_t i = 0;
    while (i < line.size()) {
        while (i < line.size() && std::isspace(static_cast<unsigned char>(line[i]))) ++i;
        size_t j = i;
        while (j < line.size() && !std::isspace(static_cast<unsigned char>(line[j]))) ++j;
        if (j > i) fields.emplace_back(line.substr(i, j - i));
        i = j;
    }
    return fields;
}

static int64_t parse_i64_field(std::string_view field, std::string_view name) {
    std::string copy(field);
    size_t idx = 0;
    long long parsed = std::stoll(copy, &idx);
    if (idx != copy.size()) {
        throw std::runtime_error("invalid integer in BED " + std::string(name) + ": " + copy);
    }
    return parsed;
}

static std::vector<Region> load_bed(const std::string &path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("could not open BED file: " + path);

    std::vector<Region> regions;
    std::string line;
    uint64_t line_no = 0;
    while (std::getline(in, line)) {
        ++line_no;
        if (line.empty() || starts_with(line, "#") || starts_with(line, "track") ||
            starts_with(line, "browser")) {
            continue;
        }

        std::vector<std::string_view> fields = split_ws(line);
        if (fields.size() < 3) continue;
        if (fields[0] == "chrom" && fields[1] == "x1") continue;

        Region region;
        region.chrom = std::string(fields[0]);
        region.start = parse_i64_field(fields[1], "start");
        region.end = parse_i64_field(fields[2], "end");
        if (fields.size() >= 6 && !fields[5].empty()) {
            region.strand = fields[5][0];
            if (region.strand != '+' && region.strand != '-' && region.strand != '.') {
                std::cerr << "Treating invalid BED strand as unstranded at line " << line_no
                          << "\n";
                region.strand = '.';
            }
        }
        if (region.end <= region.start) {
            std::cerr << "Skipping non-positive BED interval at line " << line_no << "\n";
            continue;
        }
        regions.push_back(std::move(region));
    }

    if (regions.empty()) throw std::runtime_error("BED file had no usable intervals");
    return regions;
}

static std::string chrom_alias(const std::string &chrom) {
    if (starts_with(chrom, "chr")) return chrom.substr(3);
    return "chr" + chrom;
}

static bwOverlappingIntervals_t *query_intervals(bigWigFile_t *bw, const std::string &chrom,
                                                 int64_t start, int64_t end) {
    uint32_t q_start = static_cast<uint32_t>(std::max<int64_t>(0, start));
    uint32_t q_end = static_cast<uint32_t>(std::max<int64_t>(0, end));
    if (q_end <= q_start) return nullptr;

    bwOverlappingIntervals_t *intervals =
        bwGetOverlappingIntervals(bw, chrom.c_str(), q_start, q_end);
    if (intervals && intervals->l > 0) return intervals;
    if (intervals) bwDestroyOverlappingIntervals(intervals);

    std::string alias = chrom_alias(chrom);
    return bwGetOverlappingIntervals(bw, alias.c_str(), q_start, q_end);
}

static void write_npy_f64_1d(const std::string &path, const std::vector<double> &values) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("could not create output file: " + path);

    std::string header = "{'descr': '<f8', 'fortran_order': False, 'shape': (" +
                         std::to_string(values.size()) + ",), }";
    while ((10 + header.size() + 1) % 16 != 0) header.push_back(' ');
    header.push_back('\n');

    if (header.size() > std::numeric_limits<uint16_t>::max()) {
        throw std::runtime_error("npy header too large");
    }
    uint16_t header_len = static_cast<uint16_t>(header.size());

    out.write("\x93NUMPY", 6);
    char version[2] = {1, 0};
    out.write(version, 2);
    char len_bytes[2] = {static_cast<char>(header_len & 0xff),
                         static_cast<char>((header_len >> 8) & 0xff)};
    out.write(len_bytes, 2);
    out.write(header.data(), static_cast<std::streamsize>(header.size()));
    out.write(reinterpret_cast<const char *>(values.data()),
              static_cast<std::streamsize>(values.size() * sizeof(double)));
}

static void accumulate_java_mode(const bwOverlappingIntervals_t *intervals, int64_t g_start,
                                 int64_t resolution, std::vector<double> &accumulation) {
    int64_t width = static_cast<int64_t>(accumulation.size());
    for (uint32_t i = 0; i < intervals->l; ++i) {
        double value = intervals->value[i];
        if (value <= 0) continue;
        int64_t relative = (static_cast<int64_t>(intervals->start[i]) - g_start) / resolution;
        if (relative >= 0 && relative < width) accumulation[relative] += value;
    }
}

static void accumulate_coverage_mode(const bwOverlappingIntervals_t *intervals, int64_t g_start,
                                     int64_t g_end, int64_t resolution,
                                     std::vector<double> &accumulation) {
    int64_t width = static_cast<int64_t>(accumulation.size());
    for (uint32_t i = 0; i < intervals->l; ++i) {
        double value = intervals->value[i];
        if (value <= 0) continue;

        int64_t s = std::max<int64_t>(intervals->start[i], g_start);
        int64_t e = std::min<int64_t>(intervals->end[i], g_end);
        if (e <= s) continue;

        int64_t first_bin = (s - g_start) / resolution;
        int64_t last_bin = (e - 1 - g_start) / resolution;
        for (int64_t bin = first_bin; bin <= last_bin; ++bin) {
            if (bin < 0 || bin >= width) continue;
            int64_t bin_start = g_start + bin * resolution;
            int64_t bin_end = bin_start + resolution;
            int64_t overlap = std::max<int64_t>(0, std::min(e, bin_end) - std::max(s, bin_start));
            accumulation[bin] += value * static_cast<double>(overlap) /
                                 static_cast<double>(resolution);
        }
    }
}

int main(int argc, char **argv) {
    try {
        Args args = parse_args(argc, argv);
        std::vector<Region> regions = load_bed(args.bed);

        if (bwInit(1 << 17) != 0) throw std::runtime_error("bwInit failed");
        bigWigFile_t *bw = bwOpen(args.bigwig.c_str(), nullptr, "r");
        if (!bw) throw std::runtime_error("could not open BigWig file: " + args.bigwig);

        int64_t width = 2 * args.window + 1;
        int64_t g_window = args.window * args.resolution;
        int64_t g_width = width * args.resolution;

        std::vector<double> total(static_cast<size_t>(width), 0.0);
        uint64_t contributing_regions = 0;

        for (const Region &region : regions) {
            int64_t midpoint = (region.start + region.end) / 2;
            // Place the anchor at the center of the middle bin. For even
            // resolutions, the middle bin is [midpoint-res/2, midpoint+res/2).
            int64_t g_start = midpoint - g_window - args.resolution / 2;
            int64_t g_end = g_start + g_width;

            bwOverlappingIntervals_t *intervals =
                query_intervals(bw, region.chrom, g_start, g_end);
            if (!intervals || intervals->l == 0) {
                if (intervals) bwDestroyOverlappingIntervals(intervals);
                continue;
            }

            std::vector<double> region_profile(static_cast<size_t>(width), 0.0);
            if (args.coverage_mode) {
                accumulate_coverage_mode(intervals, g_start, g_end, args.resolution,
                                         region_profile);
            } else {
                accumulate_java_mode(intervals, g_start, args.resolution, region_profile);
            }
            bwDestroyOverlappingIntervals(intervals);

            if (region.strand == '-') {
                std::reverse(region_profile.begin(), region_profile.end());
            }
            for (size_t i = 0; i < total.size(); ++i) total[i] += region_profile[i];
            ++contributing_regions;
        }

        bwClose(bw);
        bwCleanup();

        if (!args.normalize_all && contributing_regions == 0) {
            throw std::runtime_error("no BED intervals overlapped BigWig data");
        }
        uint64_t normalization_regions =
            args.normalize_all ? static_cast<uint64_t>(regions.size()) : contributing_regions;
        for (double &v : total) v /= static_cast<double>(normalization_regions);

        std::string out_path = args.out_prefix + ".npy";
        write_npy_f64_1d(out_path, total);
        std::cerr << "Wrote " << out_path << " with " << width << " bins normalized by "
                  << normalization_regions << (args.normalize_all ? " total" : " contributing")
                  << " regions (BigWig overlap=" << contributing_regions << ")\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}
