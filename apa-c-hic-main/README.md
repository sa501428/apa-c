# APA C `.hic`

Direct `.hic` APA implementation in C++ using the bundled `straw` reader.

## What it does

- Reads contact data directly from a `.hic` file at a single resolution
- Accepts multiple BEDPE loop files in one run, similar to `APA3`
- Supports aggregate normalization with `--ag-norm`
- Uses a sparse masked contact map per chromosome pair so only bins touched by loop windows are stored
- Writes either `.npy` float32 matrices or plain-text matrices, depending on the output filename

## Build

```bash
cd /Users/muhammad/Desktop/apaC/apa-c-hic-main
./build.sh
```

## Usage

```bash
./apa_hic [--ag-norm] [-k NORM] [--window W] [--min-dist BINS] [--max-dist BINS] \
  [--include-inter] [-r RESOLUTION] \
  input.hic loops1.bedpe output1.npy [loops2.bedpe output2.txt ...]
```

## Notes

- `--min-dist` and `--max-dist` are in bins, matching the `APA3` style.
- `--ag-norm` uses raw observed counts plus VC normalization vectors from the `.hic` file.
- Without `--ag-norm`, supported norms are `NONE`, `VC`, `VC_SQRT`, and `KR`.
- BEDPE parsing uses the first 6 columns only.
