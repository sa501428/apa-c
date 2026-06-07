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


@dataclass(frozen=True)
class OutputPair:
    label: str
    resolution: int
    loop_file: Path
    ours: Path
    reference: Path


@dataclass(frozen=True)
class CommandTiming:
    stage: str
    tool: str
    resolution: int | None
    elapsed_seconds: float
    reused_existing: bool
    command: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run apa_hic and clear-tools apa3 (--ag-norm --local-vc), then compare "
            "the resulting matrices with metrics and plots."
        )
    )
    parser.add_argument("--hic", required=True, type=Path, help="Input .hic file")
    parser.add_argument(
        "--loop-files",
        nargs="+",
        required=True,
        type=Path,
        help="One or more BEDPE loop files",
    )
    parser.add_argument(
        "--resolutions",
        nargs="+",
        required=True,
        type=int,
        help="One or more resolutions to validate",
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
        default=Path("/Users/muhammad/Desktop/apaC/apa-c-hic-main/apa_hic"),
        help="Path to the apa_hic binary",
    )
    parser.add_argument(
        "--clear-tools-jar",
        type=Path,
        default=Path("/Users/muhammad/Desktop/apaC/clear-tools.0.180.0.jar"),
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
    parser.add_argument("--window", type=int, default=None, help="APA window size")
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


def unique_loop_labels(loop_files: Iterable[Path]) -> list[str]:
    counts: dict[str, int] = {}
    labels: list[str] = []
    for loop_file in loop_files:
        stem = loop_file.stem
        counts[stem] = counts.get(stem, 0) + 1
        if counts[stem] == 1:
            labels.append(stem)
        else:
            labels.append(f"{stem}_{counts[stem]}")
    return labels


def build_pairs(args: argparse.Namespace) -> list[OutputPair]:
    labels = unique_loop_labels(args.loop_files)
    pairs: list[OutputPair] = []
    for resolution in args.resolutions:
        for label, loop_file in zip(labels, args.loop_files):
            pairs.append(
                OutputPair(
                    label=label,
                    resolution=resolution,
                    loop_file=loop_file,
                    ours=args.output_dir / f"apa_hic_{resolution}_{label}.npy",
                    reference=args.output_dir / f"clear_localvc_{resolution}_{label}.npy",
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
    pairs_by_resolution: dict[int, list[OutputPair]] = {}
    timings: list[CommandTiming] = []
    for pair in pairs:
        pairs_by_resolution.setdefault(pair.resolution, []).append(pair)

    for resolution, resolution_pairs in sorted(pairs_by_resolution.items()):
        if args.reuse_existing and all(pair.ours.exists() for pair in resolution_pairs):
            print(f"Reusing existing apa_hic outputs for resolution {resolution}")
            timings.append(
                CommandTiming(
                    stage="run",
                    tool="apa_hic",
                    resolution=resolution,
                    elapsed_seconds=0.0,
                    reused_existing=True,
                    command="",
                )
            )
            continue

        command = [str(args.apa_binary)]
        append_apa_hic_args(command, args, resolution)
        command.append(str(args.hic))
        for pair in resolution_pairs:
            command.extend([str(pair.loop_file), str(pair.ours)])
        elapsed_seconds = run_command(command)
        timings.append(
            CommandTiming(
                stage="run",
                tool="apa_hic",
                resolution=resolution,
                elapsed_seconds=elapsed_seconds,
                reused_existing=False,
                command=" ".join(str(piece) for piece in command),
            )
        )

    return timings


def run_clear_tools(args: argparse.Namespace, pairs: list[OutputPair]) -> list[CommandTiming]:
    pairs_by_resolution: dict[int, list[OutputPair]] = {}
    timings: list[CommandTiming] = []
    for pair in pairs:
        pairs_by_resolution.setdefault(pair.resolution, []).append(pair)

    for resolution, resolution_pairs in sorted(pairs_by_resolution.items()):
        if args.reuse_existing and all(pair.reference.exists() for pair in resolution_pairs):
            print(f"Reusing existing clear-tools outputs for resolution {resolution}")
            timings.append(
                CommandTiming(
                    stage="run",
                    tool="clear_tools_apa3",
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
        command.append(str(args.hic))
        for pair in resolution_pairs:
            command.extend([str(pair.loop_file), str(pair.reference)])
        elapsed_seconds = run_command(command)
        timings.append(
            CommandTiming(
                stage="run",
                tool="clear_tools_apa3",
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

    return {
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
        "center_ours": float(ours[ours.shape[0] // 2, ours.shape[1] // 2]),
        "center_reference": float(reference[reference.shape[0] // 2, reference.shape[1] // 2]),
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
            handle.write(f"{row['resolution']}_{row['label']}\n")
            for key, value in row.items():
                if key in {"label", "resolution"}:
                    continue
                handle.write(f"  {key}: {value}\n")
            handle.write("\n")

    manifest = {
        "hic": str(args.hic),
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
        plot_pair(pair, args.output_dir)


def write_timing_summaries(
    output_dir: Path,
    timings: list[CommandTiming],
    total_elapsed_seconds: float,
) -> None:
    rows = [
        {
            "stage": timing.stage,
            "tool": timing.tool,
            "resolution": "" if timing.resolution is None else timing.resolution,
            "elapsed_seconds": timing.elapsed_seconds,
            "reused_existing": timing.reused_existing,
            "command": timing.command,
        }
        for timing in timings
    ]

    with open(output_dir / "timing_summary.csv", "w", newline="") as handle:
        fieldnames = ["stage", "tool", "resolution", "elapsed_seconds", "reused_existing", "command"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    with open(output_dir / "timing_summary.txt", "w") as handle:
        for row in rows:
            resolution_text = row["resolution"] if row["resolution"] != "" else "all"
            handle.write(
                f"{row['tool']} [{row['stage']}] resolution={resolution_text} "
                f"elapsed_seconds={row['elapsed_seconds']:.3f} reused_existing={row['reused_existing']}\n"
            )
            if row["command"]:
                handle.write(f"  command: {row['command']}\n")
        handle.write(f"\nTotal wall clock seconds: {total_elapsed_seconds:.3f}\n")


def main() -> None:
    overall_start = time.perf_counter()
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / ".cache").mkdir(parents=True, exist_ok=True)
    (args.output_dir / ".mplconfig").mkdir(parents=True, exist_ok=True)

    pairs = build_pairs(args)
    timings = run_apa_hic(args, pairs)
    timings.extend(run_clear_tools(args, pairs))
    summary_start = time.perf_counter()
    write_summaries(args, pairs)
    timings.append(
        CommandTiming(
            stage="postprocess",
            tool="summary_and_plots",
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
            resolution=None,
            elapsed_seconds=total_elapsed_seconds,
            reused_existing=False,
            command="",
        )
    )
    write_timing_summaries(args.output_dir, timings, total_elapsed_seconds)


if __name__ == "__main__":
    main()
