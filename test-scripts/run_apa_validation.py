#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_APA_BINARY = SCRIPT_DIR / "apa-c-hic-main" / "apa_hic"
DEFAULT_CLEAR_TOOLS_JAR = Path("clear-tools.0.179.0.jar")


@dataclass(frozen=True)
class OutputPair:
    hic_label: str
    label: str
    resolution: int
    loop_file: Path
    ours: Path
    reference: Path


@dataclass(frozen=True)
class CommandTiming:
    stage: str
    tool: str
    hic_label: str
    resolution: int | None
    elapsed_seconds: float
    reused_existing: bool
    command: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run apa_hic and clear-tools apa3 (--ag-norm --local-vc) for multiple "
            ".hic files and BEDPE files, then compare matrices with metrics and plots."
        )
    )
    parser.add_argument(
        "--hic",
        nargs="+",
        required=True,
        type=Path,
        help="One or more input .hic files",
    )
    parser.add_argument(
        "--loop-files",
        nargs="+",
        required=True,
        type=Path,
        help="One or more BEDPE loop files; all are passed in one multi-output run",
    )
    parser.add_argument(
        "--resolutions",
        nargs="+",
        required=False,
        type=int,
        default=[1000, 500, 100],
        help="One or more resolutions to validate (default: 1000 500 100)",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory for generated matrices, plots, and summaries",
    )
    parser.add_argument(
        "--apa-binary",
        type=Path,
        default=DEFAULT_APA_BINARY,
        help="Path to the apa_hic binary",
    )
    parser.add_argument(
        "--clear-tools-jar",
        type=Path,
        default=DEFAULT_CLEAR_TOOLS_JAR,
        help="Path to clear-tools jar",
    )
    parser.add_argument(
        "--java-xmx",
        default="12G",
        help="Heap size passed to java as -Xmx",
    )
    parser.add_argument(
        "--apa-threads",
        type=int,
        default=None,
        help="Thread count passed to apa_hic via --threads",
    )
    parser.add_argument(
        "--clear-threads",
        type=int,
        default=None,
        help="Thread count passed to clear-tools via --threads",
    )
    parser.add_argument("--window", type=int, default=20, help="APA window size (default: 20)")
    parser.add_argument("--min-dist", type=int, default=None, help="Minimum loop distance in bins")
    parser.add_argument("--max-dist", type=int, default=None, help="Maximum loop distance in bins")
    parser.add_argument(
        "--include-inter",
        action="store_true",
        help="Include inter-chromosomal loops",
    )
    parser.add_argument(
        "--reuse-existing",
        action="store_true",
        help="Skip rerunning a tool when all expected outputs already exist",
    )
    return parser.parse_args()


def require_existing_paths(args: argparse.Namespace) -> None:
    required_paths = [
        ("apa_hic binary", args.apa_binary),
        ("clear-tools jar", args.clear_tools_jar),
        *(("input .hic", path) for path in args.hic),
        *(("BEDPE loop file", path) for path in args.loop_files),
    ]
    missing = [f"  {label}: {path}" for label, path in required_paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("Required validation inputs were not found:\n" + "\n".join(missing))


def unique_labels(paths: Iterable[Path]) -> list[str]:
    counts: dict[str, int] = {}
    labels: list[str] = []
    for path in paths:
        stem = path.stem
        counts[stem] = counts.get(stem, 0) + 1
        if counts[stem] == 1:
            labels.append(stem)
        else:
            labels.append(f"{stem}_{counts[stem]}")
    return labels


def build_pairs(args: argparse.Namespace) -> list[OutputPair]:
    hic_labels = unique_labels(args.hic)
    loop_labels = unique_labels(args.loop_files)
    pairs: list[OutputPair] = []
    for hic_label, hic_file in zip(hic_labels, args.hic):
        hic_output_dir = args.output_dir / hic_label
        for resolution in args.resolutions:
            for label, loop_file in zip(loop_labels, args.loop_files):
                pairs.append(
                    OutputPair(
                        hic_label=hic_label,
                        label=label,
                        resolution=resolution,
                        loop_file=loop_file,
                        ours=hic_output_dir / f"apa_hic_{resolution}_{label}.npy",
                        reference=hic_output_dir / f"clear_localvc_{resolution}_{label}.npy",
                    )
                )
    return pairs


def append_common_apa_args(command: list[str], args: argparse.Namespace, resolution: int) -> None:
    command.extend(["-r", str(resolution)])
    if args.window is not None:
        command.extend(["--window", str(args.window)])
    if args.min_dist is not None:
        command.extend(["--min-dist", str(args.min_dist)])
    if args.max_dist is not None:
        command.extend(["--max-dist", str(args.max_dist)])
    if args.include_inter:
        command.append("--include-inter")


def append_apa_hic_args(command: list[str], args: argparse.Namespace, resolution: int) -> None:
    append_common_apa_args(command, args, resolution)
    if args.apa_threads is not None:
        command.extend(["--threads", str(args.apa_threads)])


def append_clear_tools_args(command: list[str], args: argparse.Namespace, resolution: int) -> None:
    append_common_apa_args(command, args, resolution)
    if args.clear_threads is not None:
        command.extend(["--threads", str(args.clear_threads)])


def run_command(command: list[str]) -> float:
    command_text = " ".join(str(piece) for piece in command)
    print("Running:", command_text, flush=True)
    start = time.perf_counter()
    subprocess.run(command, check=True)
    elapsed_seconds = time.perf_counter() - start
    print(f"Completed in {elapsed_seconds:.3f}s", flush=True)
    return elapsed_seconds


def run_apa_hic(args: argparse.Namespace, pairs: list[OutputPair]) -> list[CommandTiming]:
    pairs_by_hic_resolution: dict[tuple[Path, int], list[OutputPair]] = {}
    timings: list[CommandTiming] = []
    hic_by_label = dict(zip(unique_labels(args.hic), args.hic))
    for pair in pairs:
        key = (hic_by_label[pair.hic_label], pair.resolution)
        pairs_by_hic_resolution.setdefault(key, []).append(pair)

    for (hic, resolution), resolution_pairs in sorted(pairs_by_hic_resolution.items(), key=lambda item: (str(item[0][0]), item[0][1])):
        if args.reuse_existing and all(pair.ours.exists() for pair in resolution_pairs):
            print(f"Reusing existing apa_hic outputs for resolution {resolution}")
            timings.append(
                CommandTiming(
                    stage="run",
                    tool="apa_hic",
                    hic_label=resolution_pairs[0].hic_label,
                    resolution=resolution,
                    elapsed_seconds=0.0,
                    reused_existing=True,
                    command="",
                )
            )
            continue

        command = [str(args.apa_binary)]
        append_apa_hic_args(command, args, resolution)
        command.append(str(hic))
        for pair in resolution_pairs:
            command.extend([str(pair.loop_file), str(pair.ours)])
        elapsed_seconds = run_command(command)
        timings.append(
            CommandTiming(
                stage="run",
                tool="apa_hic",
                hic_label=resolution_pairs[0].hic_label,
                resolution=resolution,
                elapsed_seconds=elapsed_seconds,
                reused_existing=False,
                command=" ".join(str(piece) for piece in command),
            )
        )

    return timings


def run_clear_tools(args: argparse.Namespace, pairs: list[OutputPair]) -> list[CommandTiming]:
    pairs_by_hic_resolution: dict[tuple[Path, int], list[OutputPair]] = {}
    timings: list[CommandTiming] = []
    hic_by_label = dict(zip(unique_labels(args.hic), args.hic))
    for pair in pairs:
        key = (hic_by_label[pair.hic_label], pair.resolution)
        pairs_by_hic_resolution.setdefault(key, []).append(pair)

    for (hic, resolution), resolution_pairs in sorted(pairs_by_hic_resolution.items(), key=lambda item: (str(item[0][0]), item[0][1])):
        if args.reuse_existing and all(pair.reference.exists() for pair in resolution_pairs):
            print(f"Reusing existing clear-tools outputs for resolution {resolution}")
            timings.append(
                CommandTiming(
                    stage="run",
                    tool="clear_tools_apa3",
                    hic_label=resolution_pairs[0].hic_label,
                    resolution=resolution,
                    elapsed_seconds=0.0,
                    reused_existing=True,
                    command="",
                )
            )
            continue

        command = [
            "java",
            f"-Xmx{args.java_xmx}",
            "-jar",
            str(args.clear_tools_jar),
            "apa3",
            "--ag-norm",
            "--local-vc",
        ]
        append_clear_tools_args(command, args, resolution)
        command.append(str(hic))
        for pair in resolution_pairs:
            command.extend([str(pair.loop_file), str(pair.reference)])
        elapsed_seconds = run_command(command)
        timings.append(
            CommandTiming(
                stage="run",
                tool="clear_tools_apa3",
                hic_label=resolution_pairs[0].hic_label,
                resolution=resolution,
                elapsed_seconds=elapsed_seconds,
                reused_existing=False,
                command=" ".join(str(piece) for piece in command),
            )
        )

    return timings


def clear_tools_style_apa_score(matrix: np.ndarray, resolution: int) -> float:
    n = matrix.shape[0]
    if resolution < 50:
        width = 50 // resolution
        center = n // 2
        numerator_value = float(
            matrix[
                center - width : center + width + 1,
                center - width : center + width + 1,
            ].mean()
        )
    else:
        numerator_value = float(matrix[n // 2, n // 2])

    denom_region = matrix[3 * n // 4 : n, 0 : n // 4]
    positive = denom_region[denom_region > 0]
    denom_value = float(positive.mean()) if positive.size else 0.0
    return (numerator_value + 1.0) / (denom_value + 1.0)


def summarize_pair(pair: OutputPair) -> dict[str, object]:
    ours = np.load(pair.ours)
    reference = np.load(pair.reference)
    diff = ours - reference
    abs_diff = np.abs(diff)
    rel_abs_diff = abs_diff / np.maximum(np.abs(reference), 1e-12)

    n = ours.shape[0]
    center = n // 2
    lower_left = ours[3 * n // 4 : n, 0 : n // 4]
    lower_left_reference = reference[3 * n // 4 : n, 0 : n // 4]
    lower_left_positive = lower_left[lower_left > 0]
    lower_left_reference_positive = lower_left_reference[lower_left_reference > 0]
    center_ours = float(ours[center, center])
    center_reference = float(reference[center, center])
    lower_left_ours = float(lower_left_positive.mean()) if lower_left_positive.size else 0.0
    lower_left_reference_value = (
        float(lower_left_reference_positive.mean())
        if lower_left_reference_positive.size
        else 0.0
    )

    return {
        "hic_label": pair.hic_label,
        "label": pair.label,
        "resolution": pair.resolution,
        "loop_file": str(pair.loop_file),
        "ours_path": str(pair.ours),
        "reference_path": str(pair.reference),
        "shape": f"{ours.shape[0]}x{ours.shape[1]}",
        "equal_exact": bool(np.array_equal(ours, reference)),
        "allclose_1e_6": bool(np.allclose(ours, reference, atol=1e-6, rtol=1e-6)),
        "sum_ours": float(ours.sum()),
        "sum_reference": float(reference.sum()),
        "center_ours": center_ours,
        "center_reference": center_reference,
        "lower_left_mean_positive_ours": lower_left_ours,
        "lower_left_mean_positive_reference": lower_left_reference_value,
        "apa_score_ours": clear_tools_style_apa_score(ours, pair.resolution),
        "apa_score_reference": clear_tools_style_apa_score(reference, pair.resolution),
        "mean_signed_diff": float(diff.mean()),
        "mean_abs_diff": float(abs_diff.mean()),
        "median_abs_diff": float(np.median(abs_diff)),
        "rmse": float(np.sqrt(np.mean(diff**2))),
        "max_abs_diff": float(abs_diff.max()),
        "mean_rel_abs_diff": float(rel_abs_diff.mean()),
        "max_rel_abs_diff": float(rel_abs_diff.max()),
    }


def plot_pair(pair: OutputPair, output_dir: Path) -> None:
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("XDG_CACHE_HOME", str(output_dir / ".cache"))
    os.environ.setdefault("MPLCONFIGDIR", str(output_dir / ".mplconfig"))

    import matplotlib.pyplot as plt

    ours = np.load(pair.ours)
    reference = np.load(pair.reference)
    diff = ours - reference

    vmax = max(float(ours.max()), float(reference.max()))
    diff_limit = max(float(np.abs(diff).max()), 1e-12)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    fig.suptitle(f"{pair.label} @ {pair.resolution} bp")

    im0 = axes[0].imshow(ours, cmap="magma", vmin=0, vmax=vmax)
    axes[0].set_title("apa_hic")
    fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

    im1 = axes[1].imshow(reference, cmap="magma", vmin=0, vmax=vmax)
    axes[1].set_title("clear-tools apa3\n(--ag-norm --local-vc)")
    fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    im2 = axes[2].imshow(diff, cmap="coolwarm", vmin=-diff_limit, vmax=diff_limit)
    axes[2].set_title("difference (ours - ref)")
    fig.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04)

    for ax in axes:
        ax.set_xlabel("col")
        ax.set_ylabel("row")

    fig.tight_layout()
    fig.savefig(output_dir / f"{pair.resolution}_{pair.label}_comparison.png", dpi=180)
    plt.close(fig)


def write_summaries(args: argparse.Namespace, pairs: list[OutputPair]) -> None:
    rows = [summarize_pair(pair) for pair in pairs]
    if not rows:
        return

    fieldnames = list(rows[0].keys())
    with open(args.output_dir / "comparison_summary.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    with open(args.output_dir / "comparison_summary.txt", "w") as handle:
        for row in rows:
            handle.write(f"{row['hic_label']}/{row['resolution']}_{row['label']}\n")
            for key, value in row.items():
                if key in {"label", "resolution"}:
                    continue
                handle.write(f"  {key}: {value}\n")
            handle.write("\n")

    manifest = {
        "hic": [str(path) for path in args.hic],
        "loop_files": [str(path) for path in args.loop_files],
        "resolutions": args.resolutions,
        "apa_binary": str(args.apa_binary),
        "clear_tools_jar": str(args.clear_tools_jar),
        "java_xmx": args.java_xmx,
        "apa_threads": args.apa_threads,
        "clear_threads": args.clear_threads,
        "window": args.window,
        "min_dist": args.min_dist,
        "max_dist": args.max_dist,
        "include_inter": args.include_inter,
        "reference_mode": "clear-tools apa3 --ag-norm --local-vc",
    }
    with open(args.output_dir / "run_manifest.json", "w") as handle:
        json.dump(manifest, handle, indent=2)

    for pair in pairs:
        plot_pair(pair, args.output_dir / pair.hic_label)


def write_timing_summaries(
    output_dir: Path,
    timings: list[CommandTiming],
    total_elapsed_seconds: float,
) -> None:
    rows = [
        {
            "stage": timing.stage,
            "tool": timing.tool,
            "hic_label": timing.hic_label,
            "resolution": "" if timing.resolution is None else timing.resolution,
            "elapsed_seconds": timing.elapsed_seconds,
            "reused_existing": timing.reused_existing,
            "command": timing.command,
        }
        for timing in timings
    ]

    with open(output_dir / "timing_summary.csv", "w", newline="") as handle:
        fieldnames = [
            "stage",
            "tool",
            "hic_label",
            "resolution",
            "elapsed_seconds",
            "reused_existing",
            "command",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    with open(output_dir / "timing_summary.txt", "w") as handle:
        for row in rows:
            resolution_text = row["resolution"] if row["resolution"] != "" else "all"
            hic_text = row["hic_label"] if row["hic_label"] else "all"
            handle.write(
                f"{row['tool']} [{row['stage']}] hic={hic_text} resolution={resolution_text} "
                f"elapsed_seconds={row['elapsed_seconds']:.3f} reused_existing={row['reused_existing']}\n"
            )
            if row["command"]:
                handle.write(f"  command: {row['command']}\n")
        handle.write(f"\nTotal wall clock seconds: {total_elapsed_seconds:.3f}\n")

    run_timings = {
        (timing.hic_label, timing.resolution, timing.tool): timing
        for timing in timings
        if timing.stage == "run" and not timing.reused_existing
    }
    comparison_rows = []
    keys = sorted(
        {
            (hic_label, resolution)
            for hic_label, resolution, _tool in run_timings
            if resolution is not None
        }
    )
    for hic_label, resolution in keys:
        ours = run_timings.get((hic_label, resolution, "apa_hic"))
        reference = run_timings.get((hic_label, resolution, "clear_tools_apa3"))
        if ours is None or reference is None:
            continue
        comparison_rows.append(
            {
                "hic_label": hic_label,
                "resolution": resolution,
                "apa_hic_seconds": ours.elapsed_seconds,
                "clear_tools_seconds": reference.elapsed_seconds,
                "apa_hic_speedup": (
                    reference.elapsed_seconds / ours.elapsed_seconds
                    if ours.elapsed_seconds > 0
                    else float("inf")
                ),
            }
        )

    if comparison_rows:
        apa_total = sum(float(row["apa_hic_seconds"]) for row in comparison_rows)
        clear_total = sum(float(row["clear_tools_seconds"]) for row in comparison_rows)
        with open(output_dir / "timing_comparison.csv", "w", newline="") as handle:
            fieldnames = [
                "hic_label",
                "resolution",
                "apa_hic_seconds",
                "clear_tools_seconds",
                "apa_hic_speedup",
            ]
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(comparison_rows)
            writer.writerow(
                {
                    "hic_label": "TOTAL",
                    "resolution": "",
                    "apa_hic_seconds": apa_total,
                    "clear_tools_seconds": clear_total,
                    "apa_hic_speedup": clear_total / apa_total if apa_total > 0 else float("inf"),
                }
            )

        with open(output_dir / "timing_comparison.txt", "w") as handle:
            handle.write(
                f"{'hic':<24} {'resolution':>10} {'apa-c (s)':>12} "
                f"{'clear (s)':>12} {'speedup':>10}\n"
            )
            for row in comparison_rows:
                handle.write(
                    f"{str(row['hic_label']):<24} {int(row['resolution']):>10} "
                    f"{float(row['apa_hic_seconds']):>12.3f} "
                    f"{float(row['clear_tools_seconds']):>12.3f} "
                    f"{float(row['apa_hic_speedup']):>9.2f}x\n"
                )
            handle.write(
                f"{'TOTAL':<24} {'':>10} {apa_total:>12.3f} "
                f"{clear_total:>12.3f} {clear_total / apa_total:>9.2f}x\n"
            )


def main() -> None:
    overall_start = time.perf_counter()
    args = parse_args()
    require_existing_paths(args)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / ".cache").mkdir(parents=True, exist_ok=True)
    (args.output_dir / ".mplconfig").mkdir(parents=True, exist_ok=True)
    for hic_label in unique_labels(args.hic):
        (args.output_dir / hic_label).mkdir(parents=True, exist_ok=True)

    pairs = build_pairs(args)
    timings = run_apa_hic(args, pairs)
    timings.extend(run_clear_tools(args, pairs))
    summary_start = time.perf_counter()
    write_summaries(args, pairs)
    timings.append(
        CommandTiming(
            stage="postprocess",
            tool="summary_and_plots",
            hic_label="",
            resolution=None,
            elapsed_seconds=time.perf_counter() - summary_start,
            reused_existing=False,
            command="",
        )
    )
    total_elapsed_seconds = time.perf_counter() - overall_start
    timings.append(
        CommandTiming(
            stage="total",
            tool="overall",
            hic_label="",
            resolution=None,
            elapsed_seconds=total_elapsed_seconds,
            reused_existing=False,
            command="",
        )
    )
    write_timing_summaries(args.output_dir, timings, total_elapsed_seconds)


if __name__ == "__main__":
    main()
