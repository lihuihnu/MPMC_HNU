"""Select affected tests; unsuccessful/cancelled runs never become checkpoints."""
import json
import os
from pathlib import Path
import re
import runpy
import subprocess
import urllib.request
import yaml
RULES = runpy.run_path('.github/ci/impact_rules.py')
MAP = json.loads(Path('.github/ci/workflow_map.json').read_text(encoding='utf-8'))

def matches(path, patterns):
    selected = False
    for pattern in patterns:
        negative = pattern.startswith('!')
        if RULES['compile_pattern'](pattern[1:] if negative else pattern).match(path):
            selected = not negative
    return selected

def select(paths, branch='main', action='synchronize'):
    result = RULES['route'](paths) if branch == 'main' else {k: False for k in RULES['ownership']}
    ad = RULES['ad_suites_for'](paths) if branch == 'main' else []
    thermo = RULES['thermo_suites_for'](paths) if result.get('thermodynamics_contracts') else []
    for path, spec in MAP['workflows'].items():
        key = spec.get('route_key')
        if not key:
            continue
        selector = spec['selector']
        if 'pull_request' not in spec['original_events']:
            result[key] = False
            continue
        eligible = matches(branch, selector.get('branches', ['**'])) and not matches(branch, selector.get('branches-ignore', []))
        # ready_for_review is an explicit cumulative revalidation of PR gates.
        types = selector.get('types', ['opened', 'reopened', 'synchronize'])
        eligible = eligible and (action in types or action == 'ready_for_review')
        result[key] = eligible and any(matches(p, selector.get('paths', ['**'])) and not matches(p, selector.get('paths-ignore', [])) for p in paths)
    if any(p.startswith('.github/ci/') for p in paths):
        print('CI governance changed: rely on governance regressions and path-specific selectors; do not fan out all science gates.')
    if result.get('sw92_profile_c_sensitivity'):
        result['sw92_thermodynamics'] = True
        result['sw92_profile_c_phase_set'] = True
    if result.get('sw92_phase_assigned_no_w'):
        result['sw92_family_vle'] = True
    return result, ad, thermo

CORE_ROUTER_JOB_ROUTES = {
    'flow-core': ['flow_core'],
    'flow-discretization': ['flow_discretization'],
    'flow-discretization-petsc': ['flow_discretization_petsc'],
    'pt-flash-backend': ['pt_flash_backend'],
    'pt-stability': ['pt_stability'],
    'pt-split': ['pt_split'],
    'pr76-pt-continuation': ['pr76_pt_continuation'],
    'pr76-max3': ['pr76_max3'],
    'cpa-stability': ['cpa_stability'],
    'cpa-split': ['cpa_split'],
    'cpa-max3': ['cpa_max3'],
    'cpa-baseline': ['cpa_baseline'],
    'cpa-pt-phase': ['cpa_pt_phase'],
    'cpa-physical-validation': ['cpa_physical_validation'],
    'thermodynamics': ['thermodynamics_contracts'],
    'physics-closure': ['physics_closure'],
    'model-configuration': ['model_configuration'],
    'sw92-thermodynamics': ['sw92_thermodynamics'],
    'sw92-family-vle': ['sw92_family_vle'],
    'sw92-profile-c-pt': ['sw92_profile_c_pt'],
    'sw92-profile-c-phase-set': ['sw92_profile_c_phase_set'],
    'sw92-profile-c-sensitivity': ['sw92_profile_c_sensitivity'],
    'sw92-physics-closure': ['sw92_physics_closure'],
    'sw92-phase-assigned-joint': ['sw92_phase_assigned_joint'],
    'sw92-phase-assigned-h-side-witness': ['sw92_phase_assigned_h_side_witness'],
    'sw92-phase-assigned-no-w': ['sw92_phase_assigned_no_w'],
    'sw92-phase-assigned-three-phase-topology': [
        'sw92_phase_assigned_three_phase',
        'sw92_phase_assigned_three_phase_closure',
        'sw92_phase_assigned_boundary',
        'sw92_profile_c_pt',
        'sw92_profile_c_phase_set',
    ],
}

def router_job_changes(before, head):
    def read(rev):
        proc = subprocess.run(
            ['git', 'show', rev + ':.github/workflows/pr_incremental_ci.yml'],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
        if proc.returncode != 0:
            return {}
        obj = yaml.safe_load(proc.stdout) or {}
        return obj.get('jobs', {})
    old, new = read(before), read(head)
    return sorted(name for name in set(old) | set(new) if old.get(name) != new.get(name))

def apply_router_job_changes(result, before, head):
    legacy = {}
    for spec in MAP['workflows'].values():
        key = spec.get('route_key')
        if not key:
            continue
        for job in spec.get('central_hashes', {}):
            legacy[job] = key
    changed = router_job_changes(before, head)
    force_ad = False
    force_thermo = False
    for job in changed:
        if job in ('impact', 'result'):
            continue
        if job in legacy:
            result[legacy[job]] = True
            continue
        if job == 'ad':
            force_ad = True
            continue
        if job == 'thermodynamics':
            result['thermodynamics_contracts'] = True
            force_thermo = True
            continue
        routes = CORE_ROUTER_JOB_ROUTES.get(job)
        if routes is None:
            raise RuntimeError('Changed central job lacks audited route mapping: ' + job)
        for route in routes:
            result[route] = True
    return changed, force_ad, force_thermo

def cmd(*args):
    return subprocess.check_output(['git', *args], text=True).strip()
def ancestor(a, b):
    return subprocess.run(['git', 'merge-base', '--is-ancestor', a, b], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0

def main():
    event = json.loads(Path(os.environ['GITHUB_EVENT_PATH']).read_text(encoding='utf-8'))
    pr = event.get('pull_request')
    trusted = not pr or pr['head']['repo']['full_name'] == os.environ['GITHUB_REPOSITORY']
    paths = []
    if pr:
        base, head = pr['base']['sha'], pr['head']['sha']
        before = cmd('merge-base', base, head)
        mode = 'cumulative PR diff (no compatible successful checkpoint)'
        if event.get('action') == 'synchronize':
            try:
                url = 'https://api.github.com/repos/' + os.environ['GITHUB_REPOSITORY'] + '/actions/workflows/pr_incremental_ci.yml/runs?event=pull_request&status=success&per_page=100'
                req = urllib.request.Request(url, headers={'Authorization': 'Bearer ' + os.environ['GH_TOKEN'], 'Accept': 'application/vnd.github+json'})
                with urllib.request.urlopen(req, timeout=20) as response:
                    runs = json.load(response)['workflow_runs']
                for run in runs:
                    if run.get('name') != 'CI - affected tests v2' or run.get('conclusion') != 'success':
                        continue
                    if run.get('head_branch') != pr['head']['ref']:
                        continue
                    if not any(p['number'] == pr['number'] and p['base']['sha'] == base for p in run.get('pull_requests', [])):
                        continue
                    candidate = run['head_sha']
                    if candidate == head or not ancestor(candidate, head):
                        continue
                    before, mode = candidate, 'since compatible successful checkpoint'
                    break
            except Exception as error:
                print('Checkpoint lookup unavailable; using cumulative diff:', type(error).__name__)
        raw = subprocess.check_output(['git', 'diff', '--name-only', '--no-renames', '-z', before, head])
        paths = [p for p in raw.decode('utf-8').split('\0') if p]
        result, ad, thermo = select(paths, pr['base']['ref'], event.get('action', ''))
        if '.github/workflows/pr_incremental_ci.yml' in paths:
            changed_jobs, force_ad, force_thermo = apply_router_job_changes(result, before, head)
            if force_ad:
                ad = list(RULES['all_ad_suites'])
            if force_thermo:
                thermo = list(RULES['all_thermo_suites'])
            print('Changed central jobs:', json.dumps(changed_jobs))
        for path in paths:
            if path.startswith(('modules/', 'tests/', 'frontend/', 'api/')) and not path.endswith('.md'):
                single, a, t = select([path], pr['base']['ref'], event.get('action', ''))
                if not any(single.values()) and not a and not t:
                    raise RuntimeError('Executable change lacks ownership: ' + path)
        print('Validation baseline:', before, '->', head, mode)
    elif os.environ.get('GITHUB_EVENT_NAME') == 'push':
        head = event['after']
        before = ''
        try:
            url = 'https://api.github.com/repos/' + os.environ['GITHUB_REPOSITORY'] + '/actions/workflows/pr_incremental_ci.yml/runs?event=push&status=success&per_page=100'
            req = urllib.request.Request(url, headers={'Authorization': 'Bearer ' + os.environ['GH_TOKEN'], 'Accept': 'application/vnd.github+json'})
            with urllib.request.urlopen(req, timeout=20) as response:
                runs = json.load(response)['workflow_runs']
            for run in runs:
                candidate = run['head_sha']
                if run.get('name') == 'CI - affected tests v2' and run.get('conclusion') == 'success' and run.get('head_branch') == event['ref'].removeprefix('refs/heads/') and candidate != head and ancestor(candidate, head):
                    before = candidate
                    break
        except Exception as error:
            print('Push checkpoint unavailable; validate matching tracked files:', type(error).__name__)
        if re.fullmatch(r'[0-9a-f]{40}', before) and before != '0' * 40 and ancestor(before, head):
            raw = subprocess.check_output(['git', 'diff', '--name-only', '--no-renames', '-z', before, head])
        else:
            raw = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', '-z', head])
        paths = [p for p in raw.decode('utf-8').split('\0') if p]
        result = {k: False for k in RULES['ownership']}
        ad, thermo = [], []
        branch = event.get('ref', '').removeprefix('refs/heads/')
        for spec in MAP['workflows'].values():
            if 'route_key' not in spec:
                continue
            selector = spec['original_events'].get('push')
            eligible = selector is not None
            selector = selector or {}
            eligible = eligible and matches(branch, selector.get('branches', ['**'])) and not matches(branch, selector.get('branches-ignore', []))
            result[spec['route_key']] = eligible and any(matches(p, selector.get('paths', ['**'])) and not matches(p, selector.get('paths-ignore', [])) for p in paths)
    else:
        result = {k: True for k in RULES['ownership']}
        result.update({s['route_key']: False for s in MAP['workflows'].values() if 'route_key' in s})
        ad, thermo = list(RULES['all_ad_suites']), list(RULES['all_thermo_suites'])
        print('Manual central dispatch preserves migrated core suites; use named legacy workflow for its manual parameters.')
    if not trusted:
        result = {k: False for k in result}
        ad, thermo = [], []
        print('Untrusted fork: private-runner jobs blocked.')
    print('Changed paths:', json.dumps(paths))
    print('Selected gates:', json.dumps([k for k, value in result.items() if value]))
    with open(os.environ['GITHUB_OUTPUT'], 'a', encoding='utf-8') as out:
        for key, value in result.items():
            out.write(f'{key}={str(value).lower()}\n')
        out.write(f'trusted={str(trusted).lower()}\n')
        out.write(f'ad_has_work={str(bool(ad)).lower()}\n')
        out.write('ad_matrix=' + json.dumps(RULES['ad_matrix'](ad), separators=(',', ':')) + '\n')
        out.write(f'thermo_has_work={str(bool(thermo)).lower()}\n')
        out.write('thermo_matrix=' + json.dumps(RULES['thermo_matrix'](thermo), separators=(',', ':')) + '\n')
if __name__ == '__main__':
    main()
