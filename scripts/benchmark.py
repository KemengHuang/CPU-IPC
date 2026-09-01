#!/usr/bin/env python3

import argparse
import csv
import statistics
import subprocess
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Benchmark the CPU-IPC headless runner")
    parser.add_argument("--exe", required=True, type=Path)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--steps", type=int, default=1)
    parser.add_argument(
        "--scene", choices=("cloth-bunny", "twisting-mat"), default="cloth-bunny"
    )
    parser.add_argument(
        "--broad-phase", choices=("spatial-hash", "lbvh"), default="lbvh"
    )
    parser.add_argument(
        "--linear-solver",
        choices=("cholmod", "suitesparse-ldl", "eigen-cg"),
        default="suitesparse-ldl",
    )
    parser.add_argument("--output", type=Path, default=Path("Output/benchmark"))
    return parser.parse_args()


def read_metrics(path):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError(f"no metrics found in {path}")
    return rows


def main():
    args = parse_args()
    if args.repeats <= 0 or args.steps <= 0:
        raise SystemExit("--repeats and --steps must be positive")

    executable = args.exe.resolve()
    samples = []
    for repeat in range(args.repeats):
        run_output = (args.output / f"run_{repeat + 1:02d}").resolve()
        command = [
            str(executable),
            "--scene",
            args.scene,
            "--steps",
            str(args.steps),
            "--broad-phase",
            args.broad_phase,
            "--linear-solver",
            args.linear_solver,
            "--output",
            str(run_output),
        ]
        completed = subprocess.run(command, check=True, text=True, capture_output=True)
        result_lines = [line for line in completed.stdout.splitlines() if line.startswith("RESULT ")]
        if len(result_lines) != 1:
            raise RuntimeError(f"unexpected runner output:\n{completed.stdout}")

        rows = read_metrics(run_output / "metrics.csv")
        samples.append(sum(float(row["step_ms"]) for row in rows))
        print(f"run {repeat + 1}: {samples[-1]:.3f} ms | {result_lines[0]}")

    print(
        "summary: "
        f"median={statistics.median(samples):.3f} ms, "
        f"min={min(samples):.3f} ms, max={max(samples):.3f} ms, "
        f"repeats={len(samples)}, steps={args.steps}"
    )


if __name__ == "__main__":
    main()
