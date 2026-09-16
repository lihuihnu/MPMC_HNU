#!/usr/bin/env python3
"""Real browser regression for hosted Classic PR through the production edge."""

from __future__ import annotations

import argparse
import http.client
import json
import mimetypes
import pathlib
import re
import ssl
import subprocess
import threading
import time
import urllib.error
import urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any

TOKEN_A = "mpmc-hosted-web-test-user-a"
TOKEN_B = "mpmc-hosted-web-test-user-b"
INTERNAL_SPOOF = "browser-spoof-must-be-overwritten"
WEBDRIVER_ELEMENT = "element-6066-11e4-a52e-4f735466cecf"


class Observed:
    def __init__(self) -> None:
        self.lock = threading.Lock()
        self.session_id: str | None = None
        self.release_status: int | None = None
        self.release_event = threading.Event()

    def note_session(self, value: str | None) -> None:
        if not value:
            return
        with self.lock:
            self.session_id = value

    def read_session(self) -> str | None:
        with self.lock:
            return self.session_id

    def note_release(self, status: int | None) -> None:
        with self.lock:
            self.release_status = status
        self.release_event.set()


def grpc_web_status(body: bytes) -> int | None:
    offset = 0
    while offset + 5 <= len(body):
        flag = body[offset]
        length = int.from_bytes(body[offset + 1 : offset + 5], "big")
        offset += 5
        if offset + length > len(body):
            break
        payload = body[offset : offset + length]
        offset += length
        if flag & 0x80:
            match = re.search(rb"(?:^|\r\n)grpc-status:\s*(\d+)", payload)
            if match:
                return int(match.group(1))
    return None


def grpc_response_status(header: str | None, body: bytes) -> int | None:
    if header is not None and header.isdecimal():
        return int(header)
    return grpc_web_status(body)


def webdriver_request(
    port: int,
    method: str,
    path: str,
    payload: dict[str, Any] | None = None,
) -> Any:
    data = None if payload is None else json.dumps(payload).encode("utf-8")
    request = urllib.request.Request(
        f"http://127.0.0.1:{port}{path}",
        data=data,
        method=method,
        headers={"content-type": "application/json; charset=utf-8"},
    )
    try:
        with urllib.request.urlopen(request, timeout=15) as response:
            decoded = json.load(response)
    except urllib.error.HTTPError as error:
        try:
            decoded = json.loads(error.read().decode("utf-8"))
        except Exception as cause:
            raise RuntimeError(
                f"WebDriver HTTP {error.code} for {path} without a valid W3C error body"
            ) from cause
    if "value" not in decoded:
        raise RuntimeError(f"invalid WebDriver response for {path}")
    value = decoded["value"]
    if isinstance(value, dict) and value.get("error"):
        raise RuntimeError(
            f"WebDriver {value.get('error')}: {value.get('message', 'unknown error')}"
        )
    return value


class Browser:
    def __init__(self, chromedriver: str, port: int = 9515) -> None:
        self.port = port
        self.process = subprocess.Popen(
            [chromedriver, f"--port={port}", "--log-level=WARNING"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        self.session: str | None = None
        self._wait_driver()
        value = webdriver_request(
            port,
            "POST",
            "/session",
            {
                "capabilities": {
                    "alwaysMatch": {
                        "browserName": "chrome",
                        "goog:chromeOptions": {
                            "args": [
                                "--headless=new",
                                "--no-sandbox",
                                "--disable-dev-shm-usage",
                                "--window-size=1440,1200",
                            ]
                        },
                    }
                }
            },
        )
        if not isinstance(value, dict) or not isinstance(value.get("sessionId"), str):
            self.close()
            raise RuntimeError("ChromeDriver did not return a session id")
        self.session = value["sessionId"]

    def _wait_driver(self) -> None:
        for _ in range(100):
            if self.process.poll() is not None:
                output = self.process.stdout.read() if self.process.stdout else ""
                raise RuntimeError(f"ChromeDriver exited early: {output[-2000:]}")
            try:
                webdriver_request(self.port, "GET", "/status")
                return
            except Exception:
                time.sleep(0.05)
        raise RuntimeError("ChromeDriver did not become ready")

    def _path(self, suffix: str) -> str:
        if self.session is None:
            raise RuntimeError("WebDriver session is closed")
        return f"/session/{self.session}{suffix}"

    def navigate(self, url: str) -> None:
        webdriver_request(self.port, "POST", self._path("/url"), {"url": url})

    def find(self, selector: str, *, xpath: bool = False) -> str:
        value = webdriver_request(
            self.port,
            "POST",
            self._path("/element"),
            {"using": "xpath" if xpath else "css selector", "value": selector},
        )
        if not isinstance(value, dict) or not isinstance(value.get(WEBDRIVER_ELEMENT), str):
            raise RuntimeError(f"element not found: {selector}")
        return value[WEBDRIVER_ELEMENT]

    def find_all(self, selector: str) -> list[str]:
        value = webdriver_request(
            self.port,
            "POST",
            self._path("/elements"),
            {"using": "css selector", "value": selector},
        )
        if not isinstance(value, list):
            raise RuntimeError(f"invalid element list: {selector}")
        return [
            item[WEBDRIVER_ELEMENT]
            for item in value
            if isinstance(item, dict) and isinstance(item.get(WEBDRIVER_ELEMENT), str)
        ]

    def wait(self, selector: str, timeout: float = 20.0) -> str:
        deadline = time.monotonic() + timeout
        last: Exception | None = None
        while time.monotonic() < deadline:
            try:
                return self.find(selector)
            except Exception as error:
                last = error
                time.sleep(0.1)
        state = self.execute(
            "return {url: location.href, title: document.title, "
            "body: (document.body?.innerText || '').slice(0, 4000), "
            "root: (document.getElementById('root')?.innerHTML || '').slice(0, 4000)};"
        )
        raise RuntimeError(
            f"timed out waiting for {selector}: {last}; page_state={state!r}"
        )

    def click(self, element: str) -> None:
        webdriver_request(
            self.port, "POST", self._path(f"/element/{element}/click"), {}
        )

    def set_value(self, element: str, value: str) -> None:
        webdriver_request(
            self.port, "POST", self._path(f"/element/{element}/clear"), {}
        )
        webdriver_request(
            self.port,
            "POST",
            self._path(f"/element/{element}/value"),
            {"text": value},
        )

    def text(self, element: str) -> str:
        value = webdriver_request(
            self.port, "GET", self._path(f"/element/{element}/text")
        )
        return value if isinstance(value, str) else ""

    def execute(self, script: str, args: list[Any] | None = None) -> Any:
        return webdriver_request(
            self.port,
            "POST",
            self._path("/execute/sync"),
            {"script": script, "args": args or []},
        )

    def close(self) -> None:
        if self.session is not None:
            try:
                webdriver_request(self.port, "DELETE", f"/session/{self.session}")
            except Exception:
                pass
            self.session = None
        if self.process.poll() is None:
            self.process.terminate()
            try:
                self.process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self.process.kill()
                self.process.wait(timeout=5)


class ProductServer:
    def __init__(
        self,
        dist: pathlib.Path,
        edge_host: str,
        edge_port: int,
        public_ca: pathlib.Path,
        client_cert: pathlib.Path,
        client_key: pathlib.Path,
    ) -> None:
        self.dist = dist
        self.edge_host = edge_host
        self.edge_port = edge_port
        self.observed = Observed()
        context = ssl.create_default_context(cafile=str(public_ca))
        context.load_cert_chain(str(client_cert), str(client_key))
        self.tls_context = context
        outer = self

        class Handler(BaseHTTPRequestHandler):
            protocol_version = "HTTP/1.1"
            server_version = "MPMCHostedWebTest/1"

            def do_GET(self) -> None:
                self._static()

            def do_POST(self) -> None:
                if self.path.startswith("/model-api/"):
                    self._proxy_model()
                else:
                    self.send_error(404)

            def _runtime_script(self) -> bytes:
                return (
                    "window.mpmcHostedWebModelOwnership={"
                    "convention:'MPMC/model/hosted-web-ownership/v1',"
                    "baseUrl:'/model-api',"
                    "identity:{getAccessToken:async()=>({accessToken:'"
                    + TOKEN_A
                    + "'})}};"
                ).encode("utf-8")

            def _static(self) -> None:
                path = self.path.split("?", 1)[0]
                if path == "/hosted-runtime.js":
                    payload = self._runtime_script()
                    self.send_response(200)
                    self.send_header("content-type", "text/javascript; charset=utf-8")
                    self.send_header("content-length", str(len(payload)))
                    self.send_header("cache-control", "no-store")
                    self.end_headers()
                    self.wfile.write(payload)
                    return

                relative = "index.html" if path in ("/", "/index.html") else path.lstrip("/")
                candidate = (outer.dist / relative).resolve()
                try:
                    candidate.relative_to(outer.dist.resolve())
                except ValueError:
                    self.send_error(403)
                    return
                if not candidate.is_file():
                    self.send_error(404)
                    return
                payload = candidate.read_bytes()
                if candidate.name == "index.html":
                    text = payload.decode("utf-8")
                    marker = '<script type="module"'
                    if marker not in text:
                        raise RuntimeError("Vite index is missing its module script")
                    text = text.replace(
                        marker,
                        '<script src="/hosted-runtime.js"></script>' + marker,
                        1,
                    )
                    payload = text.encode("utf-8")
                content_type = mimetypes.guess_type(candidate.name)[0] or "application/octet-stream"
                self.send_response(200)
                self.send_header("content-type", content_type)
                self.send_header("content-length", str(len(payload)))
                self.send_header("cache-control", "no-store")
                self.end_headers()
                self.wfile.write(payload)

            def _proxy_model(self) -> None:
                length = int(self.headers.get("content-length", "0"))
                body = self.rfile.read(length)
                upstream_headers: dict[str, str] = {
                    "host": outer.edge_host,
                    "x-mpmc-authenticated-principal": INTERNAL_SPOOF,
                }
                for name in (
                    "content-type",
                    "x-grpc-web",
                    "grpc-timeout",
                    "x-user-agent",
                    "grpc-encoding",
                    "grpc-accept-encoding",
                    "x-mpmc-model-session",
                    "authorization",
                ):
                    value = self.headers.get(name)
                    if value is not None:
                        upstream_headers[name] = value

                if "CreateModel" in self.path:
                    outer.observed.note_session(self.headers.get("x-mpmc-model-session"))

                connection = http.client.HTTPSConnection(
                    outer.edge_host,
                    outer.edge_port,
                    context=outer.tls_context,
                    timeout=130,
                )
                try:
                    connection.request("POST", self.path, body=body, headers=upstream_headers)
                    response = connection.getresponse()
                    grpc_status_header = response.getheader("grpc-status")
                    self.send_response(response.status)
                    for name, value in response.getheaders():
                        if name.lower() in {
                            "content-type",
                            "grpc-encoding",
                            "grpc-accept-encoding",
                            "grpc-status",
                            "grpc-message",
                        }:
                            self.send_header(name, value)

                    if "OpenModelSession" in self.path:
                        self.send_header("transfer-encoding", "chunked")
                        self.end_headers()
                        try:
                            while True:
                                chunk = response.read1(8192)
                                if not chunk:
                                    break
                                self.wfile.write(f"{len(chunk):X}\r\n".encode("ascii"))
                                self.wfile.write(chunk)
                                self.wfile.write(b"\r\n")
                                self.wfile.flush()
                            self.wfile.write(b"0\r\n\r\n")
                            self.wfile.flush()
                        except (BrokenPipeError, ConnectionResetError):
                            pass
                        return

                    payload = response.read()
                    self.send_header("content-length", str(len(payload)))
                    self.end_headers()
                    self.wfile.write(payload)
                    if "ReleaseModel" in self.path:
                        outer.observed.note_release(
                            grpc_response_status(grpc_status_header, payload)
                        )
                finally:
                    connection.close()

            def log_message(self, format: str, *args: object) -> None:
                return

        self.httpd = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        self.thread = threading.Thread(target=self.httpd.serve_forever, daemon=True)

    @property
    def url(self) -> str:
        return f"http://127.0.0.1:{self.httpd.server_port}/"

    def start(self) -> None:
        self.thread.start()

    def close(self) -> None:
        self.httpd.shutdown()
        self.httpd.server_close()
        self.thread.join(timeout=5)


def wait_for_session(observed: Observed, timeout: float = 10.0) -> str:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        value = observed.read_session()
        if value:
            return value
        time.sleep(0.05)
    raise RuntimeError("browser did not publish a model session id")


def post_wrong_principal(base_url: str, session_id: str) -> int | None:
    request = urllib.request.Request(
        base_url
        + "model-api/mpmc.model_configuration.v1.ModelConfigurationService/CreateModel",
        data=b"\x00\x00\x00\x00\x00",
        method="POST",
        headers={
            "content-type": "application/grpc-web+proto",
            "x-grpc-web": "1",
            "grpc-timeout": "10S",
            "authorization": f"Bearer {TOKEN_B}",
            "x-mpmc-model-session": session_id,
        },
    )
    with urllib.request.urlopen(request, timeout=15) as response:
        body = response.read()
        return grpc_response_status(response.headers.get("grpc-status"), body)


def fill_required_model(browser: Browser) -> None:
    browser.set_value(browser.find(".pr-model-name input"), "Hosted Web Classic PR regression")
    browser.execute("document.querySelector('.pr-record-details').open=true")
    record = browser.find_all(".pr-record-details input")
    values = [
        "hosted-web-classic-pr-repository-fixture",
        "hosted-web-smoke-r1",
        "modules/pt_process/src/parameter_snapshot_supplier.cpp",
        "repository-curated fixture in current source revision",
        "pr76_parameter_source_v1",
        "copied for hosted Web Classic PR integration regression",
        "test-only repository fixture; no external usage claim",
    ]
    if len(record) < len(values):
        raise RuntimeError("hosted Web record editor is incomplete")
    for element, value in zip(record, values):
        browser.set_value(element, value)

    add = browser.find("//button[normalize-space()='Add component']", xpath=True)
    for _ in range(3):
        browser.click(add)

    component_values = (
        ("methane", "Methane", "190.555", "4595000", "0"),
        ("ethane", "Ethane", "305.4", "4880000", "0.099"),
        ("propane", "Propane", "369.825", "4248000", "0.15308"),
    )
    for index, values in enumerate(component_values):
        inputs = browser.find_all(f'[data-component-index="{index}"] input')
        if len(inputs) != 6:
            raise RuntimeError(f"component {index} input shape changed: {len(inputs)}")
        browser.set_value(inputs[0], values[0])
        browser.set_value(inputs[1], values[1])
        browser.set_value(inputs[3], values[2])
        browser.set_value(inputs[4], values[3])
        browser.set_value(inputs[5], values[4])

    pairs = browser.find_all(".pr-kij-row input")
    if len(pairs) != 3:
        raise RuntimeError(f"expected three binary interaction fields, found {len(pairs)}")
    for pair in pairs:
        browser.set_value(pair, "0")

    browser.click(browser.find("//button[normalize-space()='Apply fluid model']", xpath=True))


def solve_model(browser: Browser) -> None:
    browser.wait('[data-expert-pt-solve="true"]', timeout=30)
    browser.set_value(browser.find('input[aria-label="PR pressure MPa"]'), "1")
    browser.set_value(browser.find('input[aria-label="PR temperature K"]'), "350")
    browser.set_value(browser.find('input[aria-label="PR feed methane"]'), "0.8")
    browser.set_value(browser.find('input[aria-label="PR feed ethane"]'), "0.1")
    browser.set_value(browser.find('input[aria-label="PR feed propane"]'), "0.1")
    browser.click(browser.find("//button[normalize-space()='Run PR flash']", xpath=True))
    browser.wait('[data-expert-result="accepted"]', timeout=120)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dist", type=pathlib.Path, required=True)
    parser.add_argument("--edge-host", default="api.example.test")
    parser.add_argument("--edge-port", type=int, default=8444)
    parser.add_argument("--public-ca", type=pathlib.Path, required=True)
    parser.add_argument("--client-cert", type=pathlib.Path, required=True)
    parser.add_argument("--client-key", type=pathlib.Path, required=True)
    parser.add_argument("--chromedriver", required=True)
    args = parser.parse_args()

    if not (args.dist / "index.html").is_file():
        raise RuntimeError("hosted Web Vite dist is missing index.html")

    product = ProductServer(
        args.dist.resolve(),
        args.edge_host,
        args.edge_port,
        args.public_ca,
        args.client_cert,
        args.client_key,
    )
    product.start()
    browser: Browser | None = None
    try:
        browser = Browser(args.chromedriver)
        browser.navigate(product.url)
        shell = browser.wait('[data-desktop-product-shell="true"]', timeout=20)
        browser.wait('[data-expert-workbench="ready"]', timeout=20)
        shell_text = browser.text(shell)
        if "Classic PR" not in shell_text or "Hosted calculation · authenticated model session" not in shell_text:
            raise RuntimeError("hosted Web did not render the shared authenticated Classic PR shell")
        workspace_text = browser.text(browser.find('[data-expert-workspace="pr76"]'))
        if "Hosted calculation · authenticated session" not in workspace_text:
            raise RuntimeError("hosted Web workspace did not publish authenticated-session status")

        fill_required_model(browser)
        browser.wait('[data-expert-pt-solve="true"]', timeout=30)
        session_id = wait_for_session(product.observed)
        wrong_status = post_wrong_principal(product.url, session_id)
        if wrong_status != 7:
            raise RuntimeError(
                f"cross-principal model access was not denied: grpc-status={wrong_status}"
            )

        solve_model(browser)
        browser.click(browser.find('[data-product-mode="pt"]'))
        if not product.observed.release_event.wait(15):
            raise RuntimeError("Classic PR unmount did not dispatch ReleaseModel")
        if product.observed.release_status != 0:
            raise RuntimeError(
                f"ReleaseModel did not complete successfully: {product.observed.release_status}"
            )

        print(
            "HOSTED_WEB_PRODUCT_SMOKE_OK "
            "shared_classic_pr_ui=true cross_principal_denied=true release=true",
            flush=True,
        )
        return 0
    finally:
        if browser is not None:
            browser.close()
        product.close()


if __name__ == "__main__":
    raise SystemExit(main())
