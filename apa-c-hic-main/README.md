# APA C `.hic`

Direct `.hic` APA implementation in C++ using the bundled `straw` reader.

## What it does

- Reads contact data directly from a `.hic` file at a single resolution
- Accepts multiple BEDPE loop files in one run, similar to `APA3`
- Always uses aggregate normalization
- Processes chromosome pairs in parallel with thread-local accumulators
- Uses a record-driven ROI grid so only loop windows hit by observed contacts are updated
- Recomputes the row/column coverage vectors directly from the observed contacts
- Writes either `.npy` float32 matrices or plain-text matrices, depending on the output filename

## Build

```bash
cd /Users/muhammad/Desktop/apaC/apa-c-hic-main
./build.sh
```

## Usage

```bash
./apa_hic [--window W] [--min-dist BINS] [--max-dist BINS] [--threads N] \
  [--include-inter] [-r RESOLUTION] \
  input.hic loops1.bedpe output1.npy [loops2.bedpe output2.txt ...]
```

## Notes

- `--min-dist` and `--max-dist` are in bins, matching the `APA3` style.
- `--threads` defaults to the machine's available hardware threads.
- Aggregate normalization is always used.
- Coverage vectors are recomputed locally from raw observed counts instead of relying on `.hic`-stored VC vectors.
- BEDPE parsing uses the first 6 columns only.
