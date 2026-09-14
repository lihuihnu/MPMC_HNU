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
  echo 'Usage: build_smoke_apk.sh --core-library <so> --out <dir>' >&2
  exit 2
fi

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
emulator_root="$repo_root/products/pt_android/emulator"
sdk_root="${ANDROID_SDK_ROOT:-${ANDROID_HOME:-}}"
compile_api="${MPMC_ANDROID_EMULATOR_API:-35}"
min_api="${MPMC_ANDROID_MIN_API:-26}"
build_tools_version="${MPMC_ANDROID_BUILD_TOOLS_VERSION:-35.0.0}"

if [[ -z "$sdk_root" ]]; then
  echo 'Android SDK root is unavailable.' >&2
  exit 1
fi
core_library="$(realpath "$core_library")"
out_dir="$(mkdir -p "$out_dir" && cd "$out_dir" && pwd)"
build_tools="$sdk_root/build-tools/$build_tools_version"
android_jar="$sdk_root/platforms/android-$compile_api/android.jar"

for tool in aapt2 d8 zipalign apksigner; do
  if [[ ! -x "$build_tools/$tool" ]]; then
    echo "Android build tool is unavailable: $build_tools/$tool" >&2
    exit 1
  fi
done
if [[ ! -f "$android_jar" || ! -f "$core_library" ]]; then
  echo 'Android platform jar or native core library is missing.' >&2
  exit 1
fi

classes_dir="$out_dir/classes"
dex_dir="$out_dir/dex"
payload_dir="$out_dir/payload"
rm -rf "$classes_dir" "$dex_dir" "$payload_dir"
mkdir -p "$classes_dir" "$dex_dir" "$payload_dir/lib/x86_64"

mapfile -t java_sources < <(find "$emulator_root/java" -type f -name '*.java' -print | sort)
if [[ ${#java_sources[@]} -ne 2 ]]; then
  echo "Expected exactly two smoke Java sources, found ${#java_sources[@]}" >&2
  exit 1
fi

javac \
  -encoding UTF-8 \
  -source 8 \
  -target 8 \
  -classpath "$android_jar" \
  -d "$classes_dir" \
  "${java_sources[@]}"

classes_jar="$out_dir/classes.jar"
jar --create --file "$classes_jar" -C "$classes_dir" .
"$build_tools/d8" \
  --min-api "$min_api" \
  --lib "$android_jar" \
  --output "$dex_dir" \
  "$classes_jar"

base_apk="$out_dir/base-unsigned.apk"
"$build_tools/aapt2" link \
  --manifest "$emulator_root/AndroidManifest.xml" \
  -I "$android_jar" \
  --min-sdk-version "$min_api" \
  --target-sdk-version "$compile_api" \
  --version-code 1 \
  --version-name '0.1.0-emulator-smoke' \
  -o "$base_apk"

cp "$dex_dir/classes.dex" "$payload_dir/classes.dex"
cp "$core_library" "$payload_dir/lib/x86_64/libmpmc_pt_android_core.so"
unsigned_apk="$out_dir/mpmc-pt-android-smoke-unsigned.apk"
cp "$base_apk" "$unsigned_apk"
(
  cd "$payload_dir"
  zip -q -r "$unsigned_apk" classes.dex lib
)

aligned_apk="$out_dir/mpmc-pt-android-smoke-aligned.apk"
"$build_tools/zipalign" -f 4 "$unsigned_apk" "$aligned_apk"

keystore="$out_dir/debug.keystore"
keytool -genkeypair -noprompt \
  -keystore "$keystore" \
  -storepass android \
  -keypass android \
  -alias androiddebugkey \
  -dname 'CN=Android Debug,O=MPMC_HNU,C=US' \
  -keyalg RSA \
  -keysize 2048 \
  -validity 10000 >/dev/null 2>&1

signed_apk="$out_dir/mpmc-pt-android-emulator-smoke.apk"
"$build_tools/apksigner" sign \
  --ks "$keystore" \
  --ks-key-alias androiddebugkey \
  --ks-pass pass:android \
  --key-pass pass:android \
  --out "$signed_apk" \
  "$aligned_apk"
"$build_tools/apksigner" verify --verbose "$signed_apk"

echo "ANDROID_SMOKE_APK_OK path=$signed_apk"
