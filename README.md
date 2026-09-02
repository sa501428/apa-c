# APA C

Direct aggregate peak analysis from `.hic`, 4DN/DCIC `.pairs`, or BAM contacts.
The `.hic` reader uses the sibling `straw` checkout; BAM support uses htslib.

## What it does

- Reads contact data directly from a `.hic`, `.pairs`, `.pairs.gz`, or `.bam` file at a single resolution
- Accepts multiple BEDPE loop files in one run, similar to `APA3`
- Uses locally calculated aggregate normalization by default, or a normalization stored in a `.hic` file on request
- Processes `.hic` chromosome pairs in parallel and streams pairs/BAM contacts once
- Uses a record-driven ROI grid so only loop windows hit by observed contacts are updated
- Recomputes the row/column VC coverage vectors directly from the observed contacts
- Writes either `.npy` float32 matrices or plain-text matrices, depending on the output filename

## Build

```bash
./apa-c-hic-main/build.sh
```

By default the build uses `../straw/C++` relative to this repository. Set
`STRAW_CPP_DIR=/path/to/straw/C++` to use another checkout.

If htslib is visible through `pkg-config`, BAM support is enabled automatically.
Without htslib the binary still supports `.hic` and pairs input and reports a clear
error if BAM is requested.

## Usage

```bash
./apa_hic [--window W] [--min-dist BINS] [--max-dist BINS] [--threads N] \
  [--include-inter] [-r RESOLUTION] [-k HIC_NORM] \
  [--input-format auto|hic|pairs|bam] \
  input.{hic,pairs,pairs.gz,bam} loops1.bedpe output1.npy \
  [loops2.bedpe output2.txt ...]
```

## Notes

- `--min-dist` and `--max-dist` are in bins, matching the `APA3` style.
- `--threads` defaults to the machine's available hardware threads. Work is
  dynamically balanced across chromosome pairs for `.hic`; for BAM it also
  controls htslib decompression threads.
- In the default aggregate-normalization path, a worker reuses one `.hic`
  stream across all data blocks in a chromosome pair.
- Pairs and BAM inputs are scanned once for all requested loop files and
  chromosome pairs. Each input contact contributes a raw count of one, and the
  same pass constructs the per-pair VC coverage vectors used by aggregate normalization.
- Standard 4DN/DCIC pairs coordinates are interpreted as 1-based. `#columns:`
  may put `chrom1`, `pos1`, `chrom2`, and `pos2` in any order. `#chromsize:`
  headers are used when present; otherwise chromosome lengths are inferred by
  a preliminary scan. Plain and gzip-compressed pairs files are supported.
- BAM contacts use the zero-based alignment starts from primary, mapped,
  non-QC-fail, non-duplicate read1 records and the corresponding mate starts.
  This works with queryname- or coordinate-sorted paired-end BAMs.
- Aggregate normalization remains the default. Its coverage vectors are recomputed locally from raw observed counts.
- Pass `-k KR`, `-k SCALE`, or `--norm <TYPE>` to accumulate contacts using
  that normalization from the `.hic` file and skip aggregate normalization.
  The requested normalization must exist for the selected resolution and
  chromosome pairs.
- `--ag-norm` explicitly selects the default aggregate-normalization behavior.
- Stored `.hic` normalizations are intentionally rejected for pairs and BAM inputs.
- `--input-format` is only needed when the filename extension is ambiguous.
- BEDPE parsing uses the first 6 columns only.

## Contact-input smoke test

After building, run:

```bash
python3 test-scripts/test_contact_inputs.py
```

The test creates equivalent pairs, gzipped pairs, and BAM fixtures and checks
that all three produce the same expected aggregate-normalized matrix. BAM is
skipped when `samtools` is not installed.
