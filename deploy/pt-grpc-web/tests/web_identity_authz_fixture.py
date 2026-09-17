#!/usr/bin/env python3
"""Synthetic external-authorizer fixture for the hosted Web product regression."""

from __future__ import annotations

import argparse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

TOKENS = {
    "Bearer mpmc-hosted-web-test-user-a": "test-user-a",
    "Bearer mpmc-hosted-web-test-user-b": "test-user-b",
}


class Handler(BaseHTTPRequestHandler):
    server_version = "MPMCTestAuthz/1"

    def do_GET(self) -> None:
        if self.path == "/healthz":
            self.send_response(200)
            self.send_header("content-type", "text/plain")
            self.end_headers()
            self.wfile.write(b"ok\n")
            return
        self._authorize()

    def do_POST(self) -> None:
        self._authorize()

    def _authorize(self) -> None:
        principal = TOKENS.get(self.headers.get("authorization", ""))
        if principal is None:
            self.send_response(401)
            self.send_header("www-authenticate", 'Bearer realm="mpmc-hosted-web-test"')
            self.send_header("cache-control", "no-store")
            self.end_headers()
            return
        self.send_response(200)
        self.send_header("x-mpmc-authenticated-principal", principal)
        self.send_header("x-envoy-auth-headers-to-remove", "authorization")
        self.send_header("cache-control", "no-store")
        self.end_headers()

    def log_message(self, format: str, *args: object) -> None:
        # Never log request headers: the fixture's tokens are test-only bearer material.
        return


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=10003)
    args = parser.parse_args()
    if not 1 <= args.port <= 65535:
        parser.error("--port must be in [1, 65535]")
    server = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    print(f"HOSTED_WEB_AUTHZ_READY port={args.port}", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
