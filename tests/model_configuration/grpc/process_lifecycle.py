"""Real native host: abrupt RPC client death and desktop control-pipe EOF."""
import concurrent.futures
import json
import secrets
import subprocess
import sys

host_exe, client_exe = sys.argv[1:]
token = secrets.token_urlsafe(32)
pool = concurrent.futures.ThreadPoolExecutor(max_workers=2)
processes = []


def line(pipe):
    value = pool.submit(pipe.readline).result(timeout=15)
    if not value:
        raise RuntimeError("child exited before handshake")
    return value.strip()


try:
    host = subprocess.Popen([host_exe, "--desktop-session-token-stdin", "--enable-model-sessions"],
                            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    processes.append(host)
    host.stdin.write(token + "\n")
    host.stdin.flush()
    ready = json.loads(line(host.stdout))
    assert ready["event"] == "pt_process_ready"
    address = "127.0.0.1:" + str(ready["selected_port"])
    holder = subprocess.Popen([client_exe, address, "hold"], stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    processes.append(holder)
    holder.stdin.write(token + "\n")
    holder.stdin.flush()
    sid, handle = line(holder.stdout), line(holder.stdout)
    holder.kill()  # No TryCancel, ReleaseModel or normal client shutdown.
    holder.wait(timeout=10)
    probe = subprocess.run([client_exe, address, "stale"], input="\n".join([token, sid, handle, ""]),
                           capture_output=True, text=True, timeout=15)
    if probe.returncode:
        raise RuntimeError("abrupt client disconnect did not reclaim its session")
    # A fresh session still creates and solves on the same host.
    next_client = subprocess.Popen([client_exe, address, "hold"], stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    processes.append(next_client)
    next_client.stdin.write(token + "\n")
    next_client.stdin.flush()
    next_sid, next_handle = line(next_client.stdout), line(next_client.stdout)
    assert next_sid != sid and next_handle != handle
    host.stdin.close()  # Actual parent control-pipe EOF stops the native host.
    host.wait(timeout=15)
    tail, diagnostic = host.stdout.read(), host.stderr.read()
    assert host.returncode == 0 and '"event":"pt_process_stopped"' in tail
    for secret in (token, sid, handle, next_sid, next_handle):
        assert secret not in tail + diagnostic
    print("PASS abrupt client disconnect, new session recovery and desktop EOF shutdown")
finally:
    for child in reversed(processes):
        if child.poll() is None:
            child.kill()
        child.wait(timeout=10)
    pool.shutdown(wait=True, cancel_futures=True)
