#!/usr/bin/env bash
set -euo pipefail

core_library=""
out_dir=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --core-library)
      core_library="$2"
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

if [[ -z "$core_library" || -z "$out_dir" ]]; then
  echo 'Usage: build_product_shell.sh --core-library <so> --out <dir>' >&2
  exit 2
fi

shell_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
core_library="$(realpath "$core_library")"
out_dir="$(mkdir -p "$out_dir" && cd "$out_dir" && pwd)"
android_root="$shell_root/android"

if [[ ! -f "$core_library" ]]; then
  echo "Android native core is missing: $core_library" >&2
  exit 1
fi
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

jni_dir="$android_root/app/src/main/jniLibs/x86_64"
mkdir -p "$jni_dir"
cp "$core_library" "$jni_dir/libmpmc_pt_android_core.so"

(
  cd "$android_root"
  ./gradlew --no-daemon assembleDebug
)

apk="$android_root/app/build/outputs/apk/debug/app-debug.apk"
if [[ ! -f "$apk" ]]; then
  echo "Capacitor debug APK is missing: $apk" >&2
  exit 1
fi
cp "$apk" "$out_dir/mpmc-pt-android-product-shell-v1-debug.apk"

echo "ANDROID_PRODUCT_SHELL_APK_OK path=$out_dir/mpmc-pt-android-product-shell-v1-debug.apk"
