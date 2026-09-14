#!/usr/bin/env bash
set -euo pipefail

apk=""
out_dir=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --apk)
      apk="$2"
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

if [[ -z "$apk" || -z "$out_dir" ]]; then
  echo 'Usage: run_emulator_smoke.sh --apk <apk> --out <dir>' >&2
  exit 2
fi

sdk_root="${ANDROID_SDK_ROOT:-${ANDROID_HOME:-}}"
api="${MPMC_ANDROID_EMULATOR_API:-35}"
if [[ -z "$sdk_root" ]]; then
  echo 'Android SDK root is unavailable.' >&2
  exit 1
fi

adb="$sdk_root/platform-tools/adb"
emulator="$sdk_root/emulator/emulator"
avdmanager="$sdk_root/cmdline-tools/latest/bin/avdmanager"
system_image="system-images;android-$api;google_apis;x86_64"
avd_name="mpmc-pt-api${api}-x86_64"
apk="$(realpath "$apk")"
out_dir="$(mkdir -p "$out_dir" && cd "$out_dir" && pwd)"
emulator_log="$out_dir/emulator.txt"
logcat_file="$out_dir/logcat.txt"
export ANDROID_AVD_HOME="$out_dir/avd"
mkdir -p "$ANDROID_AVD_HOME"

for tool in "$adb" "$emulator" "$avdmanager"; do
  if [[ ! -x "$tool" ]]; then
    echo "Android emulator tool is unavailable: $tool" >&2
    exit 1
  fi
done
if [[ ! -f "$apk" ]]; then
  echo "Smoke APK is missing: $apk" >&2
  exit 1
fi

rm -rf "$ANDROID_AVD_HOME/$avd_name.avd" "$ANDROID_AVD_HOME/$avd_name.ini"
echo no | "$avdmanager" create avd \
  --force \
  --name "$avd_name" \
  --package "$system_image" \
  --device pixel_2 >/dev/null

if ! "$emulator" -list-avds | grep -Fx "$avd_name" >/dev/null; then
  echo "Created AVD is not discoverable: $avd_name" >&2
  find "$ANDROID_AVD_HOME" -maxdepth 2 -print >&2 || true
  exit 1
fi

acceleration=(-accel off)
if [[ -e /dev/kvm ]]; then
  sudo chmod 666 /dev/kvm
  acceleration=(-accel on)
fi

"$emulator" \
  -avd "$avd_name" \
  -no-window \
  -no-audio \
  -no-boot-anim \
  -no-snapshot \
  -no-metrics \
  -gpu swiftshader_indirect \
  -camera-back none \
  -camera-front none \
  -memory 2048 \
  "${acceleration[@]}" \
  >"$emulator_log" 2>&1 &
emulator_pid=$!

cleanup() {
  "$adb" emu kill >/dev/null 2>&1 || true
  kill "$emulator_pid" >/dev/null 2>&1 || true
}
trap cleanup EXIT

show_emulator_diagnostics() {
  echo '--- emulator diagnostics ---' >&2
  echo "ANDROID_AVD_HOME=$ANDROID_AVD_HOME" >&2
  "$emulator" -list-avds >&2 || true
  "$adb" devices -l >&2 || true
  tail -n 200 "$emulator_log" >&2 || true
}

"$adb" start-server >/dev/null
if ! timeout 180 "$adb" wait-for-device; then
  echo 'Android emulator never connected to adb within 180 seconds.' >&2
  show_emulator_diagnostics
  exit 1
fi

booted=0
for _ in $(seq 1 120); do
  if [[ "$("$adb" shell getprop sys.boot_completed 2>/dev/null | tr -d '\r')" == "1" ]]; then
    booted=1
    break
  fi
  sleep 2
done
if [[ $booted -ne 1 ]]; then
  echo 'Android emulator connected to adb but did not complete boot.' >&2
  show_emulator_diagnostics
  exit 1
fi

abi="$("$adb" shell getprop ro.product.cpu.abi | tr -d '\r')"
if [[ "$abi" != "x86_64" ]]; then
  echo "Unexpected emulator ABI: $abi" >&2
  show_emulator_diagnostics
  exit 1
fi

"$adb" install -r "$apk" >/dev/null
"$adb" shell pm path org.mpmc.ptandroid | grep -F 'package:' >/dev/null
"$adb" logcat -c
"$adb" shell am force-stop org.mpmc.ptandroid
"$adb" shell am start -W -n org.mpmc.ptandroid/.SmokeActivity >/dev/null

success=0
for _ in $(seq 1 150); do
  tagged="$({ "$adb" logcat -d -v brief -s 'MPMC_PT_ANDROID:*' '*:S'; } 2>/dev/null || true)"
  if grep -Fq 'ANDROID_JNI_PT_SERVICE_FAIL' <<<"$tagged"; then
    printf '%s\n' "$tagged" >&2
    "$adb" logcat -d -v threadtime >"$logcat_file" || true
    exit 1
  fi
  if grep -Fq 'ANDROID_JNI_PT_SERVICE_OK backends=3' <<<"$tagged"; then
    printf '%s\n' "$tagged"
    success=1
    break
  fi
  sleep 2
done

"$adb" logcat -d -v threadtime >"$logcat_file" || true
if [[ $success -ne 1 ]]; then
  echo 'Android JNI PT service smoke did not publish success before timeout.' >&2
  show_emulator_diagnostics
  exit 1
fi

"$adb" uninstall org.mpmc.ptandroid >/dev/null
if "$adb" shell pm path org.mpmc.ptandroid 2>/dev/null | grep -Fq 'package:'; then
  echo 'Smoke application remains installed after uninstall.' >&2
  exit 1
fi

echo 'ANDROID_EMULATOR_JNI_PT_SERVICE_OK abi=x86_64'
