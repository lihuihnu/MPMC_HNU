"""Export existing native fixture and exact CMake target paths for hosted Node tests."""
import json
import os
from pathlib import Path
import subprocess
import sys

manifest = Path(sys.argv[1]).resolve()
targets = json.loads(manifest.read_text(encoding="utf-8"))
host = Path(targets["host"]).resolve(strict=True)
fixture_exe = Path(targets["fixture"]).resolve(strict=True)
result = subprocess.run([str(fixture_exe), "--fixture-json"], capture_output=True,
                        text=True, check=True, timeout=30)
data = json.loads(result.stdout)
fixture = manifest.parent / "model-client-fixture.json"
fixture.write_text(json.dumps(data), encoding="utf-8")
with open(os.environ["GITHUB_ENV"], "a", encoding="utf-8") as output:
    output.write(f"MPMC_MODEL_CLIENT_HOST={host}\nMPMC_MODEL_CLIENT_FIXTURE={fixture}\n")
print("Prepared attributed native fixture and host paths for live shared-client tests")
