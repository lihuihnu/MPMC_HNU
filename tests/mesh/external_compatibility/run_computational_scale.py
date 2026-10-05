"""Run one reproducible computational-mesh scale case outside the checkout.

The producer performs the same assertions at small/cloud and large/local sizes.
Memory is the OS process high-water mark, not a sum of estimated array sizes.
"""
import argparse
import ctypes
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time


def windows_peak(process):
    from ctypes import wintypes
    class Counters(ctypes.Structure):
        _fields_ = [('cb', wintypes.DWORD), ('faults', wintypes.DWORD)] + [
            (name, ctypes.c_size_t) for name in ('peak_working_set', 'working_set', 'peak_paged',
            'paged', 'peak_nonpaged', 'nonpaged', 'pagefile', 'peak_pagefile', 'private')]
    counters = Counters()
    counters.cb = ctypes.sizeof(counters)
    query = ctypes.WinDLL('psapi', use_last_error=True).GetProcessMemoryInfo
    query.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
    query.restype = wintypes.BOOL
    if not query(wintypes.HANDLE(int(process._handle)), ctypes.byref(counters), counters.cb):
        return None
    return int(counters.peak_working_set)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--producer', type=Path, required=True)
    parser.add_argument('--format', choices=('gmsh', 'vtu', 'grdecl'), required=True)
    parser.add_argument('--extent', type=int, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--timeout', type=float, default=1800)
    args = parser.parse_args()
    producer = args.producer.resolve(strict=True)
    directory = args.output_dir.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    command = [str(producer), '--computational-scale', args.format, str(args.extent), args.format]
    report = dict(status='failed', started_utc=datetime.now(timezone.utc).isoformat(),
                  platform=platform.platform(), python=sys.version, cpu=os.environ.get('PROCESSOR_IDENTIFIER', platform.processor()),
                  logical_cpu_count=os.cpu_count(), producer_sha256=hashlib.sha256(producer.read_bytes()).hexdigest(),
                  command=command, cwd=str(directory), timeout_seconds=args.timeout)
    start = time.monotonic()
    peak = 0
    try:
        with (directory/'run.log').open('w', encoding='utf-8') as log:
            process = subprocess.Popen(command, cwd=directory, stdout=log, stderr=subprocess.STDOUT)
            while process.poll() is None:
                if os.name == 'nt':
                    peak = max(peak, windows_peak(process) or 0)
                if time.monotonic()-start > args.timeout:
                    process.kill()
                    process.wait()
                    raise TimeoutError('producer exceeded explicit timeout')
                time.sleep(.05)
            if os.name == 'nt':
                last_peak = windows_peak(process)
                peak = max(peak, last_peak or 0)
                report['memory_measurement'] = 'GetProcessMemoryInfo PeakWorkingSetSize'
                report['post_exit_query_succeeded'] = last_peak is not None
            else:
                import resource
                peak = int(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
                if sys.platform != 'darwin':
                    peak *= 1024
                report['memory_measurement'] = 'getrusage child maxrss; one producer process'
        report['exit_code'] = process.returncode
        if process.returncode != 0:
            raise RuntimeError('producer failed; see run.log')
        marker = f'[PASS] computational.mesh.scale format={args.format} cells={args.extent**3}'
        if marker not in (directory/'run.log').read_text(encoding='utf-8'):
            raise RuntimeError('producer verification marker missing')
        report['metrics'] = dict(line.split('=', 1) for line in
                                 (directory/(args.format+'.metrics.txt')).read_text().splitlines())
        report['status'] = 'passed'
    except Exception as error:
        report['error'] = str(error)
    finally:
        report['wall_seconds'] = time.monotonic()-start
        report['peak_process_memory_bytes'] = peak or None
        report['files_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in directory.iterdir() if p.is_file()}
        (directory/'result.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({key: report.get(key) for key in
                     ('status', 'wall_seconds', 'peak_process_memory_bytes', 'error')}))
    return 0 if report['status'] == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
