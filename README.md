# APA C `.hic`

Direct `.hic` APA implementation in C++ using the sibling `straw` checkout.

## What it does

- Reads contact data directly from a `.hic` file at a single resolution
- Accepts multiple BEDPE loop files in one run, similar to `APA3`
- Uses locally calculated aggregate normalization by default, or a normalization stored in the `.hic` file on request
- Processes chromosome pairs in parallel with thread-local accumulators
- Uses a record-driven ROI grid so only loop windows hit by observed contacts are updated
- Recomputes the row/column coverage vectors directly from the observed contacts
- Writes either `.npy` float32 matrices or plain-text matrices, depending on the output filename

## Build

```bash
./apa-c-hic-main/build.sh
```

By default the build uses `../straw/C++` relative to this repository. Set
`STRAW_CPP_DIR=/path/to/straw/C++` to use another checkout.

## Usage

```bash
./apa_hic [--window W] [--min-dist BINS] [--max-dist BINS] [--threads N] \
  [--include-inter] [-r RESOLUTION] [-k HIC_NORM] \
  input.hic loops1.bedpe output1.npy [loops2.bedpe output2.txt ...]
```

## Notes

- `--min-dist` and `--max-dist` are in bins, matching the `APA3` style.
- `--threads` defaults to the machine's available hardware threads. Work is
  dynamically balanced across chromosome pairs without nested thread pools.
- In the default aggregate-normalization path, a worker reuses one `.hic`
  stream across all data blocks in a chromosome pair.
- Aggregate normalization remains the default. Its coverage vectors are recomputed locally from raw observed counts.
- Pass `-k KR`, `-k SCALE`, or `--norm <TYPE>` to accumulate contacts using
  that normalization from the `.hic` file and skip aggregate normalization.
  The requested normalization must exist for the selected resolution and
  chromosome pairs.
- `--ag-norm` explicitly selects the default aggregate-normalization behavior.
- BEDPE parsing uses the first 6 columns only.
