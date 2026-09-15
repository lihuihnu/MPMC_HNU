"""Ephemeral software-test mTLS identities, generated on official CI runners."""
import pathlib
import shutil
import subprocess
import sys

root = pathlib.Path(sys.argv[1]).resolve()
root.mkdir(parents=True, exist_ok=True)
openssl = shutil.which("openssl")
if not openssl:
    raise RuntimeError("runner OpenSSL executable is required for mTLS tests")


def run(*args):
    subprocess.run([openssl, *args], cwd=root, check=True, capture_output=True)


run("req", "-x509", "-newkey", "rsa:2048", "-nodes", "-sha256", "-days", "1",
    "-subj", "/CN=model-session-test-ca", "-keyout", "ca.key", "-out", "ca.crt")
for name in ("server", "client-a", "client-b"):
    # Both clients share a CN to expose incorrect CN-only principal binding.
    cn = "localhost" if name == "server" else "same-client-name"
    run("req", "-newkey", "rsa:2048", "-nodes", "-sha256", "-subj", "/CN=" + cn,
        "-keyout", name + ".key", "-out", name + ".csr")
    extension = "subjectAltName=DNS:localhost\nextendedKeyUsage=serverAuth\n" if name == "server" else "extendedKeyUsage=clientAuth\n"
    (root / (name + ".ext")).write_text(extension)
    run("x509", "-req", "-sha256", "-days", "1", "-in", name + ".csr", "-CA", "ca.crt",
        "-CAkey", "ca.key", "-CAcreateserial", "-extfile", name + ".ext", "-out", name + ".crt")
print("Generated ephemeral session-test certificates")
