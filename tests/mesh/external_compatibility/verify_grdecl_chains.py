#!/usr/bin/env python3
"""XTGeo-generated GRDECL -> MPMC -> XTGeo, including inactive-cell properties.

Synthetic rectilinear grids only; no reservoir physics or fault/NNC claim.
XTGeo 4.26.0's pinned internal GRDECL reader additionally exposes the raw logical
COORD/ZCORN arrays. Public Grid/GridProperty APIs verify geometry and properties.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys

import numpy as np
import xtgeo
from xtgeo.grid3d._grdecl_grid import GrdeclGrid


DIMENSIONS = (2, 2, 2)
FIELDS = ("PORO", "PERMX", "PERMY", "PERMZ")


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def close(actual, expected, label):
    # Decimal roundoff only: no large absolute floor for permeability in m².
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=0, err_msg=label)


def flat(values):
    return np.asarray(values).ravel(order="F")


def write_properties(path, fields, rle=False):
    for name, values in fields.items():
        prop = xtgeo.GridProperty(ncol=2, nrow=2, nlay=2, name=name, values=values)
        prop.to_file(path, fformat="grdecl", append=True, dtype=np.float64,
                     fmt="%.17g", rle=rle)


def read_external(path):
    logical = GrdeclGrid.from_file(path)
    grid = xtgeo.grid_from_file(path, fformat="grdecl")
    actnum = flat(grid.get_actnum().values).astype(int)
    # Use a fully active geometry copy so the public property reader exposes
    # values at inactive logical cells too. No ACTNUM state is changed on disk.
    full = grid.copy()
    full.set_actnum(xtgeo.GridProperty(full, name="ACTNUM", discrete=True, values=1))
    properties = {name: flat(xtgeo.gridproperty_from_file(
        path, fformat="grdecl", name=name, grid=full).values)
        for name in FIELDS}
    return dict(dimensions=grid.dimensions, coord=logical.coord, zcorn=logical.zcorn,
                actnum=actnum, fields=properties,
                volumes=flat(full.get_bulk_volume(asmasked=False).values),
                corners=np.asarray([full.get_xyz_cell_corners((i, j, k), activeonly=False)
                                    for k in (1, 2) for j in (1, 2) for i in (1, 2)]))


def check_external(actual, expected):
    require(actual["dimensions"] == expected["dimensions"], "SPECGRID dimensions")
    for key in ("coord", "zcorn", "corners", "volumes"):
        close(actual[key], expected[key], key)
    require(np.array_equal(actual["actnum"], expected["actnum"]), "ACTNUM binding")
    for name in FIELDS:
        close(actual["fields"][name], expected["fields"][name], name + " logical binding")


def check_import(audit, expected, coordinate_scale, permeability_scale):
    require(tuple(audit["dimensions"]) == DIMENSIONS, "import dimensions")
    close(audit["coord_m"], expected["coord"] * coordinate_scale, "import COORD SI")
    close(audit["zcorn_m"], expected["zcorn"] * coordinate_scale, "import ZCORN SI")
    require(np.array_equal(audit["actnum"], expected["actnum"]), "import ACTNUM")
    active = np.flatnonzero(expected["actnum"])
    require(audit["active_logical_ids"] == (active + 1).tolist(), "active logical ID binding")
    close(audit["active_volumes_m3"], expected["volumes"][active] * coordinate_scale**3,
          "active volumes SI")
    for name in FIELDS:
        values = expected["fields"][name] * (1 if name == "PORO" else permeability_scale)
        close(audit["fields"][name], values, name + " import SI")
        close(audit["active_fields"][name], values[active], name + " active binding SI")


def negative_controls(directory, name, output, expected, properties, has_hole):
    rejected = {}
    for control in ("coord", "zcorn", *FIELDS, *(("actnum",) if has_hole else ())):
        path = directory / (name + "_bad_" + control + ".grdecl")
        logical = GrdeclGrid.from_file(output)
        fields = {key: value.copy() for key, value in properties.items()}
        if control == "coord":
            logical.coord[0] += 1
        elif control == "zcorn":
            logical.zcorn[0] += 1
        elif control == "actnum":
            # Same count, different logical cell identity.
            logical.actnum[[0, 2]] = logical.actnum[[2, 0]]
        else:
            values = flat(fields[control]).copy()
            values[[0, 1]] = values[[1, 0]]
            fields[control] = values.reshape(DIMENSIONS, order="F")
        logical.to_file(path, fileformat="grdecl")
        write_properties(path, fields)
        # Every negative remains externally readable. Only comparison with the
        # independent source binding may reject it, never an unrelated parse error.
        actual = read_external(path)
        try:
            check_external(actual, expected)
        except AssertionError as error:
            reason = str(error)
            require(control.upper() in reason.upper(), "wrong negative-control reason: " + reason)
            rejected[control] = reason
        else:
            raise AssertionError("negative control escaped: " + control)
    return rejected


def run_case(directory, producer, name, has_hole, coordinate_scale, permeability_scale):
    path = directory / (name + "_input.grdecl")
    grid = xtgeo.create_box_grid(DIMENSIONS, origin=(10, 20, 30), increment=(2, 3, 4))
    active = np.ones(8, dtype=np.int32)
    if has_hole:
        active[2] = 0
    grid.set_actnum(xtgeo.GridProperty(grid, name="ACTNUM", discrete=True,
                                     values=active.reshape(DIMENSIONS, order="F")))
    index = np.arange(8).reshape(DIMENSIONS, order="F")
    properties = {"PORO": .125 + index / 64, "PERMX": 10. + index,
                  "PERMY": 20. + 2 * index, "PERMZ": 30. + 3 * index}
    grid.to_file(path, fformat="grdecl", rle=has_hole)
    write_properties(path, properties, rle=has_hole)
    expected = read_external(path)
    # Independent analytic rectangular-cell oracle, rather than only equality
    # between two calls to the same external reader.
    close(expected["volumes"], np.full(8, 24.), "analytic cell volume")
    require(np.array_equal(expected["actnum"], active), "generator ACTNUM")
    for field, values in properties.items():
        close(expected["fields"][field], flat(values), "generator " + field)
    stem = name + "_output"
    run = subprocess.run([str(producer), "--convert-grdecl", path.name, stem,
                          str(coordinate_scale), str(permeability_scale)],
                         cwd=directory, capture_output=True, text=True, timeout=60)
    require(run.returncode == 0, "converter failed: " + run.stdout + run.stderr)
    audit = json.loads((directory / (stem + ".import.json")).read_text())
    check_import(audit, expected, coordinate_scale, permeability_scale)
    output = directory / (stem + ".grdecl")
    actual = read_external(output)
    check_external(actual, expected)
    require((directory / (stem + ".report")).read_text().splitlines() == ["lossless"],
            "GRDECL native chain must report lossless")
    rejected = negative_controls(directory, name, output, expected, properties, has_hole)
    print("[PASS] independent.grdecl." + name)
    return dict(input=path.name, output=output.name, coordinate_scale_to_m=coordinate_scale,
                permeability_scale_to_m2=permeability_scale, dimensions=list(DIMENSIONS),
                actnum=active.tolist(), logical_volumes_file_units_cubed=actual["volumes"].tolist(),
                active_logical_ids=audit["active_logical_ids"],
                active_volumes_m3=audit["active_volumes_m3"],
                fields={key: value.tolist() for key, value in actual["fields"].items()},
                negative_controls=rejected, converter_stdout=run.stdout)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--producer", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    directory = args.output_dir.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    producer = args.producer.resolve(strict=True)
    repository = Path(__file__).resolve().parents[3]
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repository,
                          capture_output=True, text=True, check=True, timeout=30).stdout.strip()
    status = subprocess.run(["git", "status", "--porcelain"], cwd=repository,
                            capture_output=True, text=True, check=True, timeout=30).stdout
    requirements = Path(__file__).with_name("requirements-grdecl-reader.txt")
    reader = Path(sys.modules[GrdeclGrid.__module__].__file__)
    report = dict(status="failed", started_utc=datetime.now(timezone.utc).isoformat(),
                  platform=platform.platform(), python=sys.version, xtgeo=xtgeo.__version__,
                  checkout_head=head, checkout_status=status,
                  requirements_sha256=hashlib.sha256(requirements.read_bytes()).hexdigest(),
                  xtgeo_grdecl_reader_sha256=hashlib.sha256(reader.read_bytes()).hexdigest(),
                  producer_sha256=hashlib.sha256(producer.read_bytes()).hexdigest(),
                  oracle_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), cases={})
    try:
        require(xtgeo.__version__ == "4.26.0", "requires pinned XTGeo 4.26.0")
        for case in (("all_active_si", False, 1., 1.), ("inactive_scaled", True, 2.5, 1e-15)):
            report["cases"][case[0]] = run_case(directory, producer, *case)
        report.update(status="passed", chains=2, conversion_reports=2, negative_controls=13)
        print("[PASS] independent.mesh.grdecl_chains chains=2 reports=2 negative_controls=13")
    except Exception as error:
        report["error"] = str(error)
        print("[FAIL] " + str(error), file=sys.stderr)
        return 1
    finally:
        report["output_sha256"] = {p.relative_to(directory).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                                   for p in sorted(directory.rglob("*")) if p.is_file()}
        (directory / "result.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
