#!/usr/bin/env bash
set -euo pipefail

core_library_x86_64=""
core_library_arm64_v8a=""
out_dir=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --core-library-x86_64)
      core_library_x86_64="$2"
      shift 2
      ;;
    --core-library-arm64-v8a)
      core_library_arm64_v8a="$2"
      shift 2
      ;;
    --out)
      out_dir="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 2
      ;;
  esac
done

if [[ -z "$core_library_x86_64" || -z "$core_library_arm64_v8a" || -z "$out_dir" ]]; then
  echo 'Usage: build_product_shell.sh --core-library-x86_64 <so> --core-library-arm64-v8a <so> --out <dir>' >&2
  exit 2
fi

shell_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
core_library_x86_64="$(realpath "$core_library_x86_64")"
core_library_arm64_v8a="$(realpath "$core_library_arm64_v8a")"
out_dir="$(mkdir -p "$out_dir" && cd "$out_dir" && pwd)"
android_root="$shell_root/android"

for entry in \
  "x86_64:$core_library_x86_64" \
  "arm64-v8a:$core_library_arm64_v8a"; do
  abi="${entry%%:*}"
  library="${entry#*:}"
  if [[ ! -f "$library" ]]; then
    echo "Android native core is missing for $abi: $library" >&2
    exit 1
  fi
done
if [[ ! -f "$shell_root/dist/index.html" ]]; then
  echo 'Android React shell dist is missing; build the shell web bundle first.' >&2
  exit 1
fi
if [[ ! -x "$shell_root/node_modules/.bin/cap" ]]; then
  echo 'Capacitor CLI is not installed in the isolated Android shell.' >&2
  exit 1
fi

rm -rf "$android_root"
(
  cd "$shell_root"
  ./node_modules/.bin/cap add android
)

variables="$android_root/variables.gradle"
if [[ ! -f "$variables" ]]; then
  echo 'Generated Capacitor Android variables.gradle is missing.' >&2
  exit 1
fi
python3 - "$variables" <<'PY'
from pathlib import Path
import re
import sys

path = Path(sys.argv[1])
text = path.read_text(encoding="utf-8")
updated, count = re.subn(r"(?m)^\s*minSdkVersion\s*=\s*\d+\s*$", "    minSdkVersion = 26", text)
if count != 1:
    raise SystemExit(f"Expected one Capacitor minSdkVersion assignment, found {count}")
path.write_text(updated, encoding="utf-8")
PY

generated_package="$android_root/app/src/main/java/org/mpmc/ptandroid"
mkdir -p "$generated_package"
cp "$shell_root/android-src/org/mpmc/ptandroid/MainActivity.java" "$generated_package/MainActivity.java"
cp "$shell_root/android-src/org/mpmc/ptandroid/MpmcPtPlugin.java" "$generated_package/MpmcPtPlugin.java"
cp "$shell_root/android-src/org/mpmc/ptandroid/NativeBridge.java" "$generated_package/NativeBridge.java"

jni_x86_64_dir="$android_root/app/src/main/jniLibs/x86_64"
jni_arm64_v8a_dir="$android_root/app/src/main/jniLibs/arm64-v8a"
mkdir -p "$jni_x86_64_dir" "$jni_arm64_v8a_dir"
cp "$core_library_x86_64" "$jni_x86_64_dir/libmpmc_pt_android_core.so"
cp "$core_library_arm64_v8a" "$jni_arm64_v8a_dir/libmpmc_pt_android_core.so"

(
  cd "$android_root"
  ./gradlew --no-daemon assembleDebug
)

apk="$android_root/app/build/outputs/apk/debug/app-debug.apk"
if [[ ! -f "$apk" ]]; then
  echo "Capacitor debug APK is missing: $apk" >&2
  exit 1
fi

python3 - "$apk" <<'PY'
from pathlib import Path
import sys
import zipfile

apk = Path(sys.argv[1])
required = {
    "lib/x86_64/libmpmc_pt_android_core.so",
    "lib/arm64-v8a/libmpmc_pt_android_core.so",
}
with zipfile.ZipFile(apk) as archive:
    names = set(archive.namelist())
missing = sorted(required - names)
if missing:
    raise SystemExit(f"Universal Product Shell APK is missing native ABI payloads: {missing}")
print("ANDROID_PRODUCT_SHELL_ABI_OK abis=arm64-v8a,x86_64")
PY

output_apk="$out_dir/mpmc-pt-android-product-shell-v1-universal-debug.apk"
cp "$apk" "$output_apk"

echo "ANDROID_PRODUCT_SHELL_APK_OK path=$output_apk"
