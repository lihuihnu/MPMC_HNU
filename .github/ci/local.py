"""Local companion to GitHub CI: cumulative planning and existing CMake/CTest suites.

This is not a GitHub workflow interpreter or a replacement for the CI matrix.
Only audited native suites are executable. Other selected gates remain explicit.
"""
import argparse
from contextlib import contextmanager, redirect_stdout
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import re
import runpy
import shutil
import subprocess
import sys
import uuid

import yaml

ROOT = Path(__file__).resolve().parents[2]
NATIVE = {
    'legacy_mesh_core': ('mesh_core.yml', 'tests/mesh/core'),
    'legacy_discretization_core': ('discretization_core.yml', 'tests/discretization/core'),
    'flow_core': ('flow_core.yml', 'tests/flow/core'),
    'flow_discretization': ('flow_discretization.yml', 'tests/flow_discretization/core'),
}
AD = {
    'arithmetic': ('.', 'mpmc_ad_tests', '^ad[.]dual$'),
    'math': ('tests/ad/math', 'mpmc_ad_math_tests', '^ad[.]math[.]'),
    'jacobian': ('tests/ad/jacobian', 'mpmc_ad_jacobian_tests', '^ad[.]jacobian[.]'),
    'runtime': ('tests/ad/runtime', 'mpmc_ad_runtime_tests', '^ad[.]runtime[.]'),
}
THERMO = {
    'contracts': 'mpmc_thermodynamics_contract_tests',
    'pr76': 'mpmc_thermodynamics_pr76_tests',
    'pr76_mixture': 'mpmc_thermodynamics_pr76_mixture_tests',
    'pr76_pt': 'mpmc_thermodynamics_pr76_pt_tests',
}
SUITES = ['governance', *NATIVE, *(f'ad.{s}' for s in AD), *(f'thermo.{s}' for s in THERMO)]


@contextmanager
def repo_directory():
    previous = Path.cwd()
    os.chdir(ROOT)
    try:
        yield
    finally:
        os.chdir(previous)


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT)


def changed_paths(base):
    """Compare merge-base to index/worktree, including both rename sides and untracked."""
    baseline = git('merge-base', base, 'HEAD').decode().strip()
    raw = git('diff', '--name-only', '--no-renames', '-z', baseline, '--')
    raw += git('ls-files', '--others', '--exclude-standard', '-z')
    return baseline, sorted({p.decode('utf-8') for p in raw.split(b'\0') if p})


def normalize_paths(paths):
    result = []
    for path in paths:
        path = path.replace('\\', '/')
        if not path or path.startswith('/') or ':' in path or '..' in path.split('/'):
            raise ValueError(f'Expected repository-relative path: {path}')
        result.append(path.removeprefix('./'))
    return sorted(set(result))


def make_plan(paths):
    paths = normalize_paths(paths)
    with repo_directory(), redirect_stdout(io.StringIO()):
        planner = runpy.run_path(str(ROOT / '.github/ci/plan.py'))
        gates, ad, thermo = planner['select'](paths, branch='main', action='ready_for_review')
        unknown = []
        for path in paths:
            if path.startswith(('modules/', 'tests/', 'frontend/', 'api/')) and not path.endswith('.md'):
                selected, a, t = planner['select']([path], branch='main', action='ready_for_review')
                if not any(selected.values()) and not a and not t:
                    unknown.append(path)
    if unknown:
        raise ValueError('Executable change lacks cloud ownership: ' + ', '.join(unknown))
    selected = sorted(k for k, value in gates.items() if value)
    local = ['governance', *(k for k in selected if k in NATIVE),
             *(f'ad.{s}' for s in ad), *(f'thermo.{s}' for s in thermo)]
    # Thermodynamics is expanded into its original suites above.
    cloud_only = [k for k in selected if k not in NATIVE and k != 'thermodynamics_contracts']
    # Cloud also compares router job semantics across revisions. Path selection
    # alone cannot reproduce that for uncommitted workflow edits; fail closed.
    if '.github/workflows/pr_incremental_ci.yml' in paths:
        cloud_only.append('central_router_semantic_validation')
    return {'paths': paths, 'selected_cloud_gates': selected, 'local_suites': local,
            'cloud_only_gates': cloud_only,
            'scope': 'Local native configuration only; other CI platforms/modes remain unverified.'}


def workflow_text(filename):
    obj = yaml.safe_load((ROOT / '.github/workflows' / filename).read_text(encoding='utf-8'))
    return '\n'.join(str(step.get('run', '')) for job in obj['jobs'].values()
                     for step in job.get('steps', []))


def native_recipe(suite):
    """Read audited target lists/CTest expressions from the existing owning workflow."""
    filename, source = NATIVE[suite]
    text = workflow_text(filename).replace('\\\n', ' ')
    targets = re.findall(r'--target\s+([a-zA-Z0-9_ \t]+?)\s+--config', text)
    patterns = re.findall(r"-R\s+['\"]([^'\"]+)['\"]", text)
    sanitizers = re.findall(r'-D(MPMC_\w+_ENABLE_SANITIZERS)=', text)
    if len(targets) != 1 or len(patterns) != 1 or len(sanitizers) != 1 or f'-S {source}' not in text:
        raise ValueError(f'Local recipe requires review after workflow change: {filename}')
    return source, targets[0].split(), patterns[0], sanitizers[0]


def commands(suite, build_root, config, jobs):
    # MSBuild adds long target/tlog names below this directory. Keep a stable,
    # collision-resistant short name; the report retains the readable suite.
    build = build_root / hashlib.sha256(suite.encode()).hexdigest()[:12]
    def cmake(source, targets, pattern, sanitizer, suffix='', extra=()):
        folder = Path(str(build) + suffix)
        return [
            ['cmake', '-S', str(ROOT / source), '-B', str(folder),
             f'-DCMAKE_BUILD_TYPE={config}', f'-D{sanitizer}=OFF', *extra],
            ['cmake', '--build', str(folder), '--config', config, '--parallel', str(jobs),
             *(['--target', *targets] if targets else [])],
            ['ctest', '--test-dir', str(folder), '-C', config, '-R', pattern,
             '--verbose', '--no-tests=error'],
        ]
    if suite == 'governance':
        return [[sys.executable, '-B', str(ROOT / '.github/ci/verify_workflows.py')]]
    if suite in NATIVE:
        return cmake(*native_recipe(suite))
    if suite.startswith('ad.'):
        kind = suite.split('.')[1]
        source, target, pattern = AD[kind]
        result = cmake(source, [target], pattern, 'MPMC_AD_ENABLE_SANITIZERS')
        if kind == 'arithmetic':
            result += cmake('tests/ad/standalone', [], '.', 'MPMC_AD_ENABLE_SANITIZERS',
                            '-standalone', (f'-DMPMC_AD_SOURCE_DIR={ROOT / "modules/ad"}',))
            result.append([sys.executable, '-B', str(ROOT / 'tests/build/root_libraries/verify.py'),
                           '--config', config, '--sanitizer', 'OFF'])
        return result
    if suite.startswith('thermo.'):
        kind = suite.split('.')[1]
        result = cmake(f'tests/thermodynamics/{kind}', [THERMO[kind]],
                       f'^thermo[.]{kind}[.]', 'MPMC_THERMO_ENABLE_SANITIZERS')
        if kind == 'pr76_pt':
            result.append([sys.executable, '-B', str(ROOT / 'tests/thermodynamics/pr76_pt/reference_decimal.py')])
        return result
    raise ValueError(f'Unsupported local suite: {suite}')


def output_root(path):
    path = Path(path).resolve()
    if path == ROOT or ROOT in path.parents:
        raise ValueError('Build/log root must be outside the repository checkout')
    return path


def snapshot():
    status = git('status', '--porcelain=v1', '-z')
    digest = hashlib.sha256(git('diff', '--binary', 'HEAD', '--'))
    for raw in sorted(git('ls-files', '--others', '--exclude-standard', '-z').split(b'\0')):
        if raw:
            digest.update(raw)
            file = ROOT / raw.decode('utf-8')
            if file.is_file():
                digest.update(file.read_bytes())
    return {'head': git('rev-parse', 'HEAD').decode().strip(),
            'dirty': bool(status), 'worktree_diff_sha256': digest.hexdigest()}


def doctor():
    result = {'python': sys.version.split()[0], 'python_path': sys.executable,
              'platform': platform.platform(), 'pyyaml': yaml.__version__, 'tools': {}}
    for tool in ('git', 'cmake', 'ctest'):
        path = shutil.which(tool)
        result['tools'][tool] = {'path': path}
        if path:
            proc = subprocess.run([path, '--version'], capture_output=True, text=True)
            result['tools'][tool]['version'] = proc.stdout.splitlines()[0] if proc.stdout else ''
    if os.name == 'nt':
        vswhere = Path(os.environ.get('ProgramFiles(x86)', 'C:/Program Files (x86)')) / 'Microsoft Visual Studio/Installer/vswhere.exe'
        result['visual_studio'] = subprocess.check_output(
            [str(vswhere), '-latest', '-products', '*', '-requires', 'Microsoft.VisualStudio.Component.VC.Tools.x86.x64',
             '-property', 'installationPath']).decode('utf-8', errors='replace').strip() if vswhere.exists() else ''
    else:
        result['cxx'] = shutil.which(os.environ.get('CXX', 'c++'))
    result['ready_for_native_configure'] = (all(v['path'] for v in result['tools'].values())
                                           and yaml.__version__ == '6.0.2'
                                           and bool(result.get('visual_studio') or result.get('cxx')))
    return result


def execute(plan, build_root, config, jobs):
    run_dir = build_root / 'runs' / (datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ-') + uuid.uuid4().hex[:8])
    run_dir.mkdir(parents=True)
    temporary = build_root / 'tmp'
    temporary.mkdir(exist_ok=True)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONUTF8='1',
               TMP=str(temporary), TEMP=str(temporary), TMPDIR=str(temporary))
    report = {'started_at': datetime.now(timezone.utc).isoformat(), 'repository': str(ROOT),
              'before': snapshot(), 'environment': doctor(), 'config': config, 'sanitizers': 'OFF',
              'plan': plan, 'results': [], 'status': 'running'}
    report_path = run_dir / 'result.json'
    def save():
        report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    save()
    try:
        for suite in plan['local_suites']:
            entry = {'suite': suite, 'status': 'running', 'commands': []}
            report['results'].append(entry)
            log = run_dir / f'{suite}.log'
            entry['log'] = str(log)
            with log.open('wb') as stream:
                for command in commands(suite, build_root / config, config, jobs):
                    print(f'RUN {suite}: {command[0]}', flush=True)
                    entry['commands'].append(command)
                    stream.write((json.dumps(command, ensure_ascii=False) + '\n').encode('utf-8'))
                    stream.flush()
                    proc = subprocess.run(command, cwd=ROOT, env=env, stdout=stream, stderr=subprocess.STDOUT)
                    if proc.returncode:
                        entry.update(status='failed', exit_code=proc.returncode)
                        raise RuntimeError(f'{suite} failed; see {log}')
            entry['status'] = 'passed'
            save()
        report['status'] = 'partial' if plan['cloud_only_gates'] else 'passed_local_scope'
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError, KeyboardInterrupt) as error:
        report.update(status='failed', error=str(error))
        for entry in report['results']:
            if entry['status'] == 'running':
                entry['status'] = 'failed'
    finally:
        report['after'] = snapshot()
        if report['after'] != report['before']:
            report.update(status='failed', error='Source changed during validation; evidence is not a stable snapshot.')
        report['finished_at'] = datetime.now(timezone.utc).isoformat()
        save()
    print(json.dumps({'status': report['status'], 'report': str(report_path)}, ensure_ascii=False))
    return 0 if report['status'] == 'passed_local_scope' else (2 if report['status'] == 'partial' else 1)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['doctor', 'plan', 'run'])
    parser.add_argument('--base', default='origin/main', help='Local Git ref; fetch it explicitly before planning')
    parser.add_argument('--path', action='append', help='Preview a repository-relative path (plan only)')
    parser.add_argument('--suite', action='append', choices=SUITES, help='Explicit local scope; not a full affected-test claim')
    parser.add_argument('--build-root', help='Required for run; outside repository')
    parser.add_argument('--config', choices=['Debug', 'Release'], default='Release')
    parser.add_argument('--jobs', type=int, default=2)
    args = parser.parse_args(argv)
    if args.jobs < 1 or args.jobs > 32:
        parser.error('--jobs must be between 1 and 32')
    if args.action == 'doctor':
        info = doctor()
        print(json.dumps(info, ensure_ascii=False, indent=2))
        return 0 if info['ready_for_native_configure'] else 1
    if args.path and args.action != 'plan':
        parser.error('--path is a preview only; run uses actual Git changes or explicit --suite')
    if args.suite:
        plan = {'local_suites': list(dict.fromkeys(args.suite)), 'cloud_only_gates': [],
                'scope': 'Explicit suites only; not complete affected validation or CI matrix parity.'}
    else:
        baseline, paths = (None, args.path) if args.path else changed_paths(args.base)
        plan = make_plan(paths)
        plan.update(base_ref=args.base, merge_base=baseline)
    if args.action == 'plan':
        print(json.dumps(plan, ensure_ascii=False, indent=2))
        return 0
    if not args.build_root:
        parser.error('run requires --build-root outside the repository')
    return execute(plan, output_root(args.build_root), args.config, args.jobs)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, OSError, subprocess.SubprocessError) as error:
        print(str(error), file=sys.stderr)
        sys.exit(1)
