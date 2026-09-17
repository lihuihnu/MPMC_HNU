#!/usr/bin/env python3
"""Extract exact full-workload CPA hot-path execution counts from gcov JSON.

This script does not infer scientific work from timing. It reads line execution
counts from an instrumented Release build of the unmodified production headers.
"""

from __future__ import annotations

import argparse
import gzip
import json
import pathlib
import subprocess
import sys
import tempfile


SPECS = {
    "phase_evaluations": (
        "modules/thermodynamics/include/mpmc/thermodynamics/cpa_phase.hpp",
        "    CpaPhaseState result;",
    ),
    "root_search_calls": (
        "modules/thermodynamics/include/mpmc/thermodynamics/cpa_pt_phase.hpp",
        "        CpaPtRootSet result;",
    ),
    "root_density_evaluations": (
        "modules/thermodynamics/include/mpmc/thermodynamics/cpa_pt_phase.hpp",
        "            ++result.evaluations;",
    ),
    "association_solve_calls": (
        "modules/thermodynamics/include/mpmc/thermodynamics/cpa_association.hpp",
        "    CpaAssociationResult result;",
    ),
    "association_iteration_sweeps": (
        "modules/thermodynamics/include/mpmc/thermodynamics/cpa_association.hpp",
        "        double residual = 0.0;",
    ),
    "stability_adapter_evaluations": (
        "modules/flash/include/mpmc/flash/cpa_stability.hpp",
        "        const auto roots = model_.roots(",
    ),
    "split_adapter_evaluations": (
        "modules/flash/include/mpmc/flash/cpa_split.hpp",
        "    const auto roots = model.roots(",
    ),
}


def unique_line(path: pathlib.Path, needle: str) -> int:
    matches = [
        index
        for index, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1)
        if line == needle
    ]
    if len(matches) != 1:
        raise RuntimeError(f"expected one match for {needle!r} in {path}, found {matches}")
    return matches[0]


def normalize(path: str) -> str:
    return path.replace("\\", "/")


def find_file(files: list[dict], relative: str) -> dict:
    normalized_relative = relative.replace("\\", "/")
    suffix = "/" + normalized_relative
    matches = []
    for item in files:
        candidate = normalize(item["file"])
        if candidate == normalized_relative or candidate.endswith(suffix):
            matches.append(item)
    if len(matches) != 1:
        raise RuntimeError(
            f"expected one gcov file matching {normalized_relative}, found "
            f"{[item.get('file') for item in matches]}"
        )
    return matches[0]


def line_count(file_record: dict, line_number: int) -> int:
    matches = [line for line in file_record.get("lines", []) if line.get("line_number") == line_number]
    if len(matches) != 1:
        raise RuntimeError(
            f"expected one gcov line {line_number} in {file_record.get('file')}, found {len(matches)}"
        )
    return int(matches[0].get("count", 0))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--object-dir", required=True)
    parser.add_argument("--benchmark-source", required=True)
    parser.add_argument("--source-root", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--label", required=True)
    args = parser.parse_args()

    object_dir = pathlib.Path(args.object_dir).resolve()
    benchmark_source = pathlib.Path(args.benchmark_source).resolve()
    source_root = pathlib.Path(args.source_root).resolve()
    output = pathlib.Path(args.output).resolve()

    with tempfile.TemporaryDirectory(prefix="cpa-gcov-") as temporary:
        temp = pathlib.Path(temporary)
        command = [
            "gcov",
            "--json-format",
            "--preserve-paths",
            "-o",
            str(object_dir),
            str(benchmark_source),
        ]
        completed = subprocess.run(command, cwd=temp, text=True, capture_output=True)
        if completed.returncode != 0:
            print(completed.stdout, file=sys.stderr)
            print(completed.stderr, file=sys.stderr)
            raise RuntimeError("gcov JSON generation failed")
        archives = list(temp.glob("*.gcov.json.gz"))
        if len(archives) != 1:
            raise RuntimeError(f"expected one gcov JSON archive, found {archives}")
        with gzip.open(archives[0], "rt", encoding="utf-8") as handle:
            document = json.load(handle)

    files = document.get("files", [])
    counts: dict[str, int] = {}
    evidence: dict[str, dict[str, object]] = {}
    for name, (relative, needle) in SPECS.items():
        source_path = source_root / relative
        line = unique_line(source_path, needle)
        file_record = find_file(files, relative)
        count = line_count(file_record, line)
        counts[name] = count
        evidence[name] = {"source": relative, "line": line, "needle": needle}

    if counts["phase_evaluations"] != counts["association_solve_calls"]:
        raise RuntimeError(
            "each CPA density-state phase evaluation should invoke exactly one association solve"
        )
    if counts["root_search_calls"] != (
        counts["stability_adapter_evaluations"] + counts["split_adapter_evaluations"]
    ):
        raise RuntimeError(
            "root-search calls no longer equal stability plus split adapter evaluations"
        )
    if counts["association_iteration_sweeps"] < counts["association_solve_calls"]:
        raise RuntimeError("association iteration sweep count is inconsistent")

    result = {
        "schema": "MPMC_HNU/CPA/performance-structural-counts/v1",
        "label": args.label,
        "counts": counts,
        "derived": {
            "density_evaluations_per_root_search": (
                counts["root_density_evaluations"] / counts["root_search_calls"]
                if counts["root_search_calls"] else 0.0
            ),
            "phase_evaluations_per_root_search": (
                counts["phase_evaluations"] / counts["root_search_calls"]
                if counts["root_search_calls"] else 0.0
            ),
            "association_sweeps_per_solve": (
                counts["association_iteration_sweeps"] / counts["association_solve_calls"]
                if counts["association_solve_calls"] else 0.0
            ),
        },
        "evidence": evidence,
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
