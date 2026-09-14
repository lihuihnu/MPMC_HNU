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
  echo 'Usage: run_product_shell_emulator.sh --apk <apk> --out <dir>' >&2
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
avd_name="mpmc-pt-shell-api${api}-x86_64"
apk="$(realpath "$apk")"
out_dir="$(mkdir -p "$out_dir" && cd "$out_dir" && pwd)"
emulator_log="$out_dir/emulator.txt"
logcat_file="$out_dir/logcat.txt"
dom_file="$out_dir/dom.json"
avd_home="$out_dir/avd"
mkdir -p "$avd_home"
export ANDROID_AVD_HOME="$avd_home"

for tool in "$adb" "$emulator" "$avdmanager"; do
  if [[ ! -x "$tool" ]]; then
    echo "Android emulator tool is unavailable: $tool" >&2
    exit 1
  fi
done
if [[ ! -f "$apk" ]]; then
  echo "Android product shell APK is missing: $apk" >&2
  exit 1
fi
if ! command -v node >/dev/null 2>&1; then
  echo 'Node.js is required for the WebView DevTools DOM audit.' >&2
  exit 1
fi

echo no | "$avdmanager" create avd \
  --force \
  --name "$avd_name" \
  --package "$system_image" \
  --device pixel_2 >/dev/null
if ! "$emulator" -list-avds | grep -Fxq "$avd_name"; then
  echo "Generated Android product-shell AVD is not discoverable: $avd_name" >&2
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
  "$adb" forward --remove tcp:9222 >/dev/null 2>&1 || true
  "$adb" emu kill >/dev/null 2>&1 || true
  kill "$emulator_pid" >/dev/null 2>&1 || true
}
trap cleanup EXIT

show_diagnostics() {
  echo '--- Android product shell emulator diagnostics ---' >&2
  "$adb" devices -l >&2 || true
  "$adb" logcat -d -v brief -s 'MPMC_PT_ANDROID:*' 'chromium:*' '*:S' >&2 || true
  tail -n 200 "$emulator_log" >&2 || true
}

"$adb" start-server >/dev/null
if ! timeout 180 "$adb" wait-for-device; then
  echo 'Android product shell emulator never connected to adb.' >&2
  show_diagnostics
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
  echo 'Android product shell emulator did not complete boot.' >&2
  show_diagnostics
  exit 1
fi

abi="$("$adb" shell getprop ro.product.cpu.abi | tr -d '\r')"
if [[ "$abi" != "x86_64" ]]; then
  echo "Unexpected Android product shell emulator ABI: $abi" >&2
  exit 1
fi

"$adb" install -r "$apk" >/dev/null
"$adb" shell pm path org.mpmc.ptandroid | grep -F 'package:' >/dev/null
"$adb" logcat -c
"$adb" shell am force-stop org.mpmc.ptandroid
"$adb" shell am start -W -n org.mpmc.ptandroid/.MainActivity >/dev/null

pr76_id='pr76.methane-ethane-propane.literature-r1'
sw92_id='sw92.carbon-dioxide-water.freshwater.literature-r1'
cpa_id='cpa.methanol-water-333.15k.cr1.literature-r1'
success=0
for _ in $(seq 1 180); do
  tagged="$({ "$adb" logcat -d -v brief -s 'MPMC_PT_ANDROID:*' 'chromium:*' '*:S'; } 2>/dev/null || true)"
  if grep -Fq 'ANDROID_PRODUCT_SHELL_WEB_FAIL' <<<"$tagged"; then
    printf '%s\n' "$tagged" >&2
    "$adb" logcat -d -v threadtime >"$logcat_file" || true
    exit 1
  fi
  if grep -Fq 'ANDROID_PRODUCT_SHELL_DISCOVERY_OK backends=3' <<<"$tagged" && \
     grep -Fq "ANDROID_PRODUCT_SHELL_SOLVE_OK backend=$pr76_id" <<<"$tagged" && \
     grep -Fq "ANDROID_PRODUCT_SHELL_SOLVE_OK backend=$sw92_id" <<<"$tagged" && \
     grep -Fq "ANDROID_PRODUCT_SHELL_SOLVE_OK backend=$cpa_id" <<<"$tagged"; then
    printf '%s\n' "$tagged"
    success=1
    break
  fi
  sleep 2
done

"$adb" logcat -d -v threadtime >"$logcat_file" || true
if [[ $success -ne 1 ]]; then
  echo 'Android Product Shell did not complete discovery and three JS->Capacitor->JNI solves.' >&2
  show_diagnostics
  exit 1
fi

# Android's accessibility hierarchy intentionally treats WebView as one native
# node, so uiautomator cannot prove which React text rendered. Product Shell v1
# is a debuggable engineering APK; MainActivity enables WebView DevTools only
# when FLAG_DEBUGGABLE is set. Audit the live DOM through that local-only CDP
# endpoint instead of weakening the UI-render requirement or using OCR.
webview_socket=""
for _ in $(seq 1 40); do
  webview_socket="$(
    "$adb" shell cat /proc/net/unix 2>/dev/null |
      tr -d '\r' |
      sed -n 's/.*@\(webview_devtools_remote[^ ]*\).*/\1/p' |
      head -n 1
  )"
  if [[ -n "$webview_socket" ]]; then
    break
  fi
  sleep 0.25
done
if [[ -z "$webview_socket" ]]; then
  echo 'Debug WebView DevTools socket was not published.' >&2
  show_diagnostics
  exit 1
fi

"$adb" forward --remove tcp:9222 >/dev/null 2>&1 || true
"$adb" forward tcp:9222 "localabstract:$webview_socket" >/dev/null

node --input-type=module - "$dom_file" <<'NODE'
import fs from 'node:fs';
const output = process.argv[2];

const sleep = (milliseconds) => new Promise((resolve) => setTimeout(resolve, milliseconds));

async function findPage() {
  const deadline = Date.now() + 10_000;
  while (Date.now() < deadline) {
    try {
      const response = await fetch('http://127.0.0.1:9222/json');
      if (response.ok) {
        const pages = await response.json();
        const page = pages.find(
          (candidate) => candidate.type === 'page' && candidate.webSocketDebuggerUrl,
        );
        if (page) return page;
      }
    } catch {
      // The WebView DevTools endpoint may need a moment after ADB forwarding.
    }
    await sleep(250);
  }
  throw new Error('No debuggable Product Shell WebView page was published.');
}

async function evaluate(page) {
  return await new Promise((resolve, reject) => {
    const socket = new WebSocket(page.webSocketDebuggerUrl);
    const timeout = setTimeout(() => {
      socket.close();
      reject(new Error('Timed out evaluating the Product Shell WebView DOM.'));
    }, 10_000);

    socket.addEventListener('open', () => {
      socket.send(
        JSON.stringify({
          id: 1,
          method: 'Runtime.evaluate',
          params: {
            expression:
              '({text: document.body?.innerText ?? "", smoke: document.documentElement.dataset.mpmcAndroidProductShellSmoke ?? "", rootChildren: document.getElementById("root")?.childElementCount ?? 0})',
            returnByValue: true,
          },
        }),
      );
    });
    socket.addEventListener('message', (event) => {
      const message = JSON.parse(String(event.data));
      if (message.id !== 1) return;
      clearTimeout(timeout);
      socket.close();
      if (message.error) {
        reject(new Error(`CDP Runtime.evaluate failed: ${JSON.stringify(message.error)}`));
        return;
      }
      const value = message.result?.result?.value;
      if (typeof value !== 'object' || value === null) {
        reject(new Error('CDP did not return the Product Shell DOM audit object.'));
        return;
      }
      resolve(value);
    });
    socket.addEventListener('error', () => {
      clearTimeout(timeout);
      reject(new Error('Unable to connect to the Product Shell WebView CDP socket.'));
    });
  });
}

const page = await findPage();
const value = await evaluate(page);
fs.writeFileSync(
  output,
  `${JSON.stringify({ pageUrl: page.url, pageTitle: page.title, ...value }, null, 2)}\n`,
  'utf8',
);
if (value.smoke !== 'ok') {
  throw new Error(`React Product Shell smoke marker is ${JSON.stringify(value.smoke)}.`);
}
if (!Number.isInteger(value.rootChildren) || value.rootChildren < 1) {
  throw new Error('React Product Shell root has no rendered child elements.');
}
if (
  typeof value.text !== 'string' ||
  !value.text.includes('MPMC_HNU') ||
  !value.text.includes('Model-neutral PT Flash')
) {
  throw new Error('Rendered Product Shell DOM is missing the shared React PT UI text.');
}
console.log('ANDROID_PRODUCT_SHELL_DOM_OK shared_react_ui=true');
NODE

"$adb" uninstall org.mpmc.ptandroid >/dev/null
if "$adb" shell pm path org.mpmc.ptandroid 2>/dev/null | grep -Fq 'package:'; then
  echo 'Android Product Shell remains installed after uninstall.' >&2
  exit 1
fi

echo 'ANDROID_PRODUCT_SHELL_V1_OK abi=x86_64 backends=3'
