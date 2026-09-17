#!/usr/bin/env python3
"""Summarize paired hosted-runner CPA performance measurements.

The report is intentionally descriptive. It does not create a performance gate
from hosted wall time and does not claim statistical significance from a small
runner sample.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import statistics


def tokens(line: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for token in line.strip().split():
        if "=" in token:
            key, value = token.split("=", 1)
            result[key] = value
    return result


def parse_timing(path: pathlib.Path, marker: str) -> dict[str, list[dict[str, str]]]:
    groups = {"main": [], "head": []}
    for line in path.read_text(encoding="utf-8").splitlines():
        if marker not in line:
            continue
        data = tokens(line)
        ref = data.get("ref")
        if ref not in groups:
            raise RuntimeError(f"missing/invalid ref in timing line: {line}")
        groups[ref].append(data)
    if not groups["main"] or not groups["head"]:
        raise RuntimeError(f"missing {marker} samples in {path}")
    return groups


def parse_structure(path: pathlib.Path) -> dict[str, int | float]:
    lines = [line for line in path.read_text(encoding="utf-8").splitlines()
             if "CPA_PERF_STRUCTURE" in line]
    if len(lines) != 1:
        raise RuntimeError(f"expected one structural summary in {path}, found {len(lines)}")
    data = tokens(lines[0])
    result: dict[str, int | float] = {}
    for key, value in data.items():
        if key == "checksum":
            result[key] = float(value)
        else:
            result[key] = int(value)
    return result


def stats(values: list[float]) -> dict[str, float]:
    median = statistics.median(values)
    deviations = [abs(value - median) for value in values]
    return {
        "median": median,
        "min": min(values),
        "max": max(values),
        "mad": statistics.median(deviations),
    }


def timing_metric(groups: dict[str, list[dict[str, str]]], numerator: str,
                  denominator: str | None = None) -> dict[str, object]:
    result: dict[str, object] = {}
    for ref in ("main", "head"):
        values: list[float] = []
        for sample in groups[ref]:
            value = float(sample[numerator])
            if denominator is not None:
                value /= float(sample[denominator])
            values.append(value)
        result[ref] = stats(values)
        result[f"{ref}_samples"] = values
    main_median = float(result["main"]["median"])  # type: ignore[index]
    head_median = float(result["head"]["median"])  # type: ignore[index]
    result["head_over_main"] = head_median / main_median
    return result


def relative_mad(summary: dict[str, float]) -> float:
    median = summary["median"]
    return summary["mad"] / median if median else 0.0


def fmt_ns(value: float) -> str:
    if value >= 1.0e9:
        return f"{value / 1.0e9:.3f} s"
    if value >= 1.0e6:
        return f"{value / 1.0e6:.3f} ms"
    if value >= 1.0e3:
        return f"{value / 1.0e3:.3f} us"
    return f"{value:.3f} ns"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--full-log", required=True)
    parser.add_argument("--micro-log", required=True)
    parser.add_argument("--main-structure", required=True)
    parser.add_argument("--head-structure", required=True)
    parser.add_argument("--main-coverage", required=True)
    parser.add_argument("--head-coverage", required=True)
    parser.add_argument("--main-sha", required=True)
    parser.add_argument("--head-sha", required=True)
    parser.add_argument("--markdown", required=True)
    parser.add_argument("--json", required=True)
    args = parser.parse_args()

    full_groups = parse_timing(pathlib.Path(args.full_log), "CPA_PERF_FULL")
    micro_groups = parse_timing(pathlib.Path(args.micro_log), "CPA_PERF_MICRO")
    main_structure = parse_structure(pathlib.Path(args.main_structure))
    head_structure = parse_structure(pathlib.Path(args.head_structure))
    main_coverage = json.loads(pathlib.Path(args.main_coverage).read_text(encoding="utf-8"))
    head_coverage = json.loads(pathlib.Path(args.head_coverage).read_text(encoding="utf-8"))

    timing = {
        "full_per_five_state_workload_ns": timing_metric(full_groups, "elapsed_ns"),
        "association_per_solve_ns": timing_metric(
            micro_groups, "association_ns", "association_units"),
        "phase_per_density_state_ns": timing_metric(micro_groups, "phase_ns", "phase_units"),
        "root_per_call_ns": timing_metric(micro_groups, "root_ns", "root_units"),
        "stability_adapter_per_call_ns": timing_metric(
            micro_groups, "stability_ns", "stability_units"),
        "split_adapter_per_call_ns": timing_metric(micro_groups, "split_ns", "split_units"),
    }

    structure_keys = sorted(set(main_structure) | set(head_structure))
    structural_changes = {
        key: {"main": main_structure.get(key), "head": head_structure.get(key)}
        for key in structure_keys
        if key != "checksum" and main_structure.get(key) != head_structure.get(key)
    }
    coverage_changes = {
        key: {"main": main_coverage["counts"].get(key),
              "head": head_coverage["counts"].get(key)}
        for key in sorted(set(main_coverage["counts"]) | set(head_coverage["counts"]))
        if main_coverage["counts"].get(key) != head_coverage["counts"].get(key)
    }

    full = timing["full_per_five_state_workload_ns"]
    ratio = float(full["head_over_main"])
    noise = max(
        relative_mad(full["main"]),  # type: ignore[arg-type]
        relative_mad(full["head"]),  # type: ignore[arg-type]
    )
    observational_margin = max(0.05, 3.0 * noise)
    if ratio > 1.0 + observational_margin:
        signal = "possible_regression_signal_repeat_before_claim"
    elif ratio < 1.0 - observational_margin:
        signal = "possible_improvement_signal_repeat_before_claim"
    else:
        signal = "no_clear_hosted_runner_regression_signal"

    head_assoc = float(timing["association_per_solve_ns"]["head"]["median"])  # type: ignore[index]
    head_phase = float(timing["phase_per_density_state_ns"]["head"]["median"])  # type: ignore[index]
    association_share_of_phase = head_assoc / head_phase if head_phase else math.nan
    phase_per_root = float(head_coverage["derived"]["phase_evaluations_per_root_search"])
    sweeps_per_assoc = float(head_coverage["derived"]["association_sweeps_per_solve"])
    if phase_per_root >= 100.0 and association_share_of_phase >= 0.5:
        bottleneck = "density-root scan multiplied by repeated association fixed-point solves"
    elif phase_per_root >= 100.0:
        bottleneck = "density-root scan / repeated density-state property evaluation"
    else:
        bottleneck = "no single structural multiplier established by this audit"

    report = {
        "schema": "MPMC_HNU/CPA/performance-baseline-report/v1",
        "main_sha": args.main_sha,
        "head_sha": args.head_sha,
        "timing": timing,
        "structure": {"main": main_structure, "head": head_structure,
                      "changed": structural_changes},
        "coverage": {"main": main_coverage, "head": head_coverage,
                     "changed": coverage_changes},
        "interpretation": {
            "hosted_timing_signal": signal,
            "full_head_over_main": ratio,
            "full_observational_margin": observational_margin,
            "association_share_of_direct_phase_time": association_share_of_phase,
            "phase_evaluations_per_root_search": phase_per_root,
            "association_sweeps_per_solve": sweeps_per_assoc,
            "dominant_nested_cost": bottleneck,
            "performance_gate": "none_baseline_only",
        },
    }
    pathlib.Path(args.json).write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines = [
        "# CPA flash hot-path baseline",
        "",
        f"- main: `{args.main_sha}`",
        f"- PR head: `{args.head_sha}`",
        "- runner timing interpretation: **" + signal + "**",
        "- performance gate: **none (baseline/audit only)**",
        "",
        "## Repeated Release timings",
        "",
        "| Metric | main median | head median | head/main |",
        "| --- | ---: | ---: | ---: |",
    ]
    labels = [
        ("full_per_five_state_workload_ns", "five-state full flash"),
        ("association_per_solve_ns", "association solve / call"),
        ("phase_per_density_state_ns", "density-state phase evaluation / call"),
        ("root_per_call_ns", "density-root search / call"),
        ("stability_adapter_per_call_ns", "stability adapter / call"),
        ("split_adapter_per_call_ns", "split adapter / call"),
    ]
    for key, label in labels:
        item = timing[key]
        main_median = float(item["main"]["median"])  # type: ignore[index]
        head_median = float(item["head"]["median"])  # type: ignore[index]
        lines.append(
            f"| {label} | {fmt_ns(main_median)} | {fmt_ns(head_median)} | "
            f"{float(item['head_over_main']):.4f} |"
        )

    lines += [
        "",
        "## Exact structural counts for one five-state workload",
        "",
        f"High-level result counters changed between main/head: **{bool(structural_changes)}**.",
        f"Gcov production hot-path counts changed between main/head: **{bool(coverage_changes)}**.",
        "",
        "| Gcov count | main | head |",
        "| --- | ---: | ---: |",
    ]
    for key in sorted(head_coverage["counts"]):
        lines.append(
            f"| `{key}` | {main_coverage['counts'][key]} | {head_coverage['counts'][key]} |"
        )
    lines += [
        "",
        "## Bottleneck interpretation",
        "",
        f"- phase evaluations per root search: **{phase_per_root:.2f}**",
        f"- association fixed-point sweeps per association solve: **{sweeps_per_assoc:.2f}**",
        f"- direct association time / direct density-state phase time: "
        f"**{association_share_of_phase:.3f}**",
        f"- dominant nested cost indicated by this baseline: **{bottleneck}**",
        "",
        "Hosted wall time is reported descriptively; it is not promoted to a hard CI performance gate.",
    ]
    pathlib.Path(args.markdown).write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
