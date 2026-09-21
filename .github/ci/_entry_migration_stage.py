"""Generate a single-entry candidate from an immutable source; never update a PR ref."""
import copy
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import urllib.request
import yaml

BASE = 'bcf2ce8dfaa8fd6f63749db1fac1db6dfb571a6f'
ROUTER = '.github/workflows/pr_incremental_ci.yml'
CATALOG = '.github/ci/workflow_map.json'
def git(*args):
    return subprocess.check_output(['git', *args], text=True).strip()
def load(text):
    obj = yaml.safe_load(text)
    if True in obj:
        obj['on'] = obj.pop(True)
    return obj
def dump(obj):
    return yaml.safe_dump(obj, sort_keys=False, allow_unicode=True, width=160)
def digest(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False).encode()).hexdigest()
def put(path, text):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(text, encoding='utf-8')

paths = git('ls-tree', '-r', '--name-only', BASE, '.github/workflows').splitlines()
originals = {p: git('show', BASE + ':' + p) + '\n' for p in paths if p.endswith(('.yml', '.yaml'))}
workflows = {p: load(s) for p, s in originals.items()}
root = copy.deepcopy(workflows[ROUTER])
old_impact = root['jobs']['impact']
route_step = next(s for s in old_impact['steps'] if s.get('id') == 'route')
rule_source = route_step['run'].split('event_name = os.environ.get("EVENT_NAME", "")')[0]
assert 'def route(' in rule_source and 'Thermodynamics selector regressions' in rule_source
put('.github/ci/impact_rules.py', rule_source)
registry = {'schema': 'single-entry-v2', 'migration_base': BASE, 'workflows': {}}
legacy_jobs = {}
for path, original in workflows.items():
    events = original.get('on') or {}
    assert isinstance(events, dict), (path, events)
    unexpected = set(events) - {'pull_request', 'workflow_call', 'workflow_dispatch'}
    assert not unexpected, ('special event needs explicit preservation', path, sorted(unexpected))
    spec = {'name': original.get('name', path), 'original_events': copy.deepcopy(events),
            'manual_inputs': copy.deepcopy((events.get('workflow_dispatch') or {}).get('inputs', {})),
            'permissions': copy.deepcopy(original.get('permissions', {})),
            'env': copy.deepcopy(original.get('env', {})),
            'defaults': copy.deepcopy(original.get('defaults', {})),
            'concurrency': copy.deepcopy(original.get('concurrency', {})),
            'jobs': {k: {'runner': v.get('runs-on'), 'strategy': v.get('strategy'),
                        'condition': v.get('if'), 'environment': v.get('environment'),
                        'permissions': v.get('permissions'), 'timeout_minutes': v.get('timeout-minutes')}
                     for k, v in original.get('jobs', {}).items()}}
    registry['workflows'][path] = spec
    if path == ROUTER:
        continue
    wf = copy.deepcopy(original)
    if 'pull_request' not in events:
        spec['execution'] = 'existing reusable/manual entry'
        spec['retained_hash'] = digest(wf)
        continue
    key = 'legacy_' + Path(path).stem
    spec['route_key'] = key
    spec['selector'] = copy.deepcopy(events.get('pull_request') or {})
    assert not (set(spec['selector']) - {'branches', 'branches-ignore', 'paths', 'paths-ignore', 'types'}), path
    wf['on'].pop('pull_request')
    if not wf['on']:
        wf['on']['workflow_dispatch'] = None
        spec['entry_repair'] = 'Former PR-only gate: add parameter-free manual entry; automatic body remains in central router.'
    assert wf['jobs'] == original['jobs'], ('manual jobs changed', path)
    if 'workflow_dispatch' in events:
        assert wf['on']['workflow_dispatch'] == events['workflow_dispatch'], ('manual input drift', path)
    spec['retained_hash'] = digest(wf)
    spec['execution'] = 'PR jobs hosted by central router; original manual/reusable semantics retained'
    put(path, dump(wf))
    jobs = original['jobs']
    removed = {k for k, j in jobs.items() if j.get('uses', '').endswith('/_incremental_path_gate.yml')}
    rename = {k: key + '__' + k for k in jobs if k not in removed}
    spec['central_jobs'] = rename
    spec['central_hashes'] = {}
    # PR event has no workflow_dispatch inputs. Do not inject manual defaults
    # into PR jobs: that could accidentally enable signing or release controls.
    def rewrite_string(s):
        for old, new in sorted(rename.items(), key=lambda kv: -len(kv[0])):
            s = re.sub(r'needs\.' + re.escape(old) + r'(?=[.\s\[])', 'needs.' + new, s)
        for old in removed:
            s = s.replace('needs.' + old + '.outputs.affected', "'true'")
        s = re.sub(r'(?<![\w.])inputs\.[A-Za-z_][A-Za-z0-9_-]*', "''", s)
        s = s.replace('github.workflow', "'" + original.get('name', path).replace("'", "''") + "'")
        return s
    def rewrite(obj):
        if isinstance(obj, str):
            return rewrite_string(obj)
        if isinstance(obj, list):
            return [rewrite(v) for v in obj]
        if isinstance(obj, dict):
            return {k: rewrite(v) for k, v in obj.items()}
        return obj
    for old, new in rename.items():
        source = jobs[old]
        job = rewrite(copy.deepcopy(source))
        needs = source.get('needs', [])
        if isinstance(needs, str):
            needs = [needs]
        job['needs'] = ['impact'] + [rename[n] for n in needs if n not in removed]
        condition = str(job.pop('if', 'true')).strip().removeprefix('${{').removesuffix('}}').strip()
        job['if'] = "github.event_name == 'pull_request' && needs.impact.outputs.trusted == 'true' && needs.impact.outputs." + key + " == 'true' && (" + condition + ')'
        job['name'] = original.get('name', key) + ' / ' + str(job.get('name', old))
        if 'uses' not in job:
            if original.get('env'):
                job['env'] = {**rewrite(original['env']), **job.get('env', {})}
            if original.get('defaults'):
                merged = copy.deepcopy(original['defaults'])
                for section, values in job.get('defaults', {}).items():
                    merged[section] = {**merged.get(section, {}), **values}
                job['defaults'] = rewrite(merged)
        if 'permissions' not in job and 'permissions' in original:
            job['permissions'] = copy.deepcopy(original['permissions'])
        if original.get('concurrency') and 'concurrency' not in job:
            concurrency = rewrite(copy.deepcopy(original['concurrency']))
            if isinstance(concurrency, dict):
                concurrency['group'] = str(concurrency['group']) + '-' + new
            else:
                concurrency = str(concurrency) + '-' + new
            job['concurrency'] = concurrency
        assert job.get('strategy') == rewrite(source.get('strategy')), ('matrix drift', path, old)
        assert job.get('steps') == rewrite(source.get('steps')), ('step drift', path, old)
        spec['central_hashes'][new] = digest(job)
        legacy_jobs[new] = job

planner = r'''"""Select affected tests; unsuccessful/cancelled runs never become checkpoints."""
import json
import os
from pathlib import Path
import re
import runpy
import subprocess
import urllib.request
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
        eligible = matches(branch, selector.get('branches', ['**'])) and not matches(branch, selector.get('branches-ignore', []))
        # ready_for_review is an explicit cumulative revalidation of PR gates.
        types = selector.get('types', ['opened', 'reopened', 'synchronize'])
        eligible = eligible and (action in types or action == 'ready_for_review')
        result[key] = eligible and any(matches(p, selector.get('paths', ['**'])) and not matches(p, selector.get('paths-ignore', [])) for p in paths)
    if any(p in ('.github/ci/impact_rules.py', '.github/ci/plan.py', '.github/ci/workflow_map.json') for p in paths):
        result = {k: True for k in result}
        ad = list(RULES['all_ad_suites'])
        thermo = list(RULES['all_thermo_suites'])
    return result, ad, thermo

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
        for path in paths:
            if path.startswith(('modules/', 'tests/', 'frontend/', 'api/')) and not path.endswith('.md'):
                single, a, t = select([path], pr['base']['ref'], event.get('action', ''))
                if not any(single.values()) and not a and not t:
                    raise RuntimeError('Executable change lacks ownership: ' + path)
        print('Validation baseline:', before, '->', head, mode)
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
'''
put('.github/ci/plan.py', planner)
root['name'] = 'CI - affected tests v2'
root['on'] = {'pull_request': {'types': ['opened', 'reopened', 'synchronize', 'ready_for_review']}, 'workflow_dispatch': None}
root['permissions'] = {'contents': 'read', 'actions': 'read'}
impact = copy.deepcopy(old_impact)
impact['permissions'] = {'contents': 'read', 'actions': 'read'}
impact['outputs']['trusted'] = '${{ steps.route.outputs.trusted }}'
for spec in registry['workflows'].values():
    if 'route_key' in spec:
        key = spec['route_key']
        impact['outputs'][key] = '${{ steps.route.outputs.' + key + ' }}'
impact['steps'] = [s for s in impact['steps'] if s.get('id') != 'route']
impact['steps'].extend([
    {'name': 'Install pinned metadata validator', 'run': 'python3 -m pip install --disable-pip-version-check PyYAML==6.0.2'},
    {'name': 'Verify workflow mapping and selection contracts', 'run': 'python3 .github/ci/verify_workflows.py'},
    {'name': 'Select affected gates from verified checkpoint', 'id': 'route', 'env': {'GH_TOKEN': '${{ github.token }}'}, 'run': 'python3 .github/ci/plan.py'}])
root['jobs']['impact'] = impact
for key, job in root['jobs'].items():
    if key == 'impact':
        continue
    old = str(job.get('if', 'true')).removeprefix('${{').removesuffix('}}').strip()
    job['if'] = "needs.impact.outputs.trusted == 'true' && (" + old + ')'
root['jobs'].update(legacy_jobs)
root['jobs']['result'] = {'name': 'Required CI result', 'if': '${{ always() }}',
    'needs': list(root['jobs']), 'runs-on': 'ubuntu-24.04', 'timeout-minutes': 5,
    'permissions': {'contents': 'read'}, 'steps': [{'name': 'Reject failures, cancellations and missing impact',
    'env': {'NEEDS_JSON': '${{ toJSON(needs) }}'}, 'shell': 'python', 'run': "import json, os\nneeds = json.loads(os.environ['NEEDS_JSON'])\nassert needs['impact']['result'] == 'success', 'impact did not complete'\nassert needs['impact']['outputs'].get('trusted') == 'true', 'fork validation requires trusted execution'\nfailed = {k:v['result'] for k,v in needs.items() if v['result'] not in ('success', 'skipped')}\nassert not failed, failed\nprint('Selected CI jobs completed; unselected or conditionally inapplicable jobs are skipped, not test passes.')\n"}]}
put(ROUTER, dump(root))
put(CATALOG, json.dumps(registry, ensure_ascii=False, indent=2) + '\n')
verify = r'''"""Executable entry/mapping guard. Update the map only with an audited CI change."""
import hashlib
import json
from pathlib import Path
import runpy
import yaml

def load(path):
    obj = yaml.safe_load(Path(path).read_text(encoding='utf-8'))
    if True in obj:
        obj['on'] = obj.pop(True)
    return obj

def digest(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False).encode()).hexdigest()

def main():
    catalog = json.loads(Path('.github/ci/workflow_map.json').read_text(encoding='utf-8'))
    router_path = '.github/workflows/pr_incremental_ci.yml'
    paths = {str(p) for p in Path('.github/workflows').glob('*.yml')}
    assert paths == set(catalog['workflows']), 'workflow added or removed without an entry mapping'
    root = load(router_path)
    auto = []
    for path, spec in catalog['workflows'].items():
        wf = load(path)
        events = wf.get('on') or {}
        if set(events) - {'workflow_call', 'workflow_dispatch'}:
            auto.append(path)
        if path != router_path:
            assert events, ('lost only entry', path)
            assert digest(wf) == spec['retained_hash'], ('workflow semantics differ from audited map', path)
            if 'workflow_dispatch' in spec['original_events']:
                assert events['workflow_dispatch'] == spec['original_events']['workflow_dispatch'], path
        for job, expected in spec.get('central_hashes', {}).items():
            assert digest(root['jobs'][job]) == expected, ('central job semantic drift', job)
    assert auto == [router_path], ('multiple automatic workflow entries', auto)
    assert 'result' in root['jobs'] and root['jobs']['result']['if'] == '${{ always() }}'
    # Existing selector regression vectors are run when importing the planner.
    planner = runpy.run_path('.github/ci/plan.py')
    select = planner['select']
    results, ad, thermo = select(['.github/AGENTS.md', 'tests/AGENTS.md'])
    assert not any(results.values()) and not ad and not thermo
    results, _, _ = select(['tests/flash/cpa_clapeyron_oracle/changed.json'])
    assert results['legacy_cpa_clapeyron_oracle']
    results, _, _ = select(['frontend/src/example.ts'])
    assert results['legacy_frontend']
    results, _, _ = select(['tests/flow_discretization/petsc/changed.cpp'])
    assert results['flow_discretization_petsc'] and not results['legacy_frontend']
    # Deleted/renamed files are represented by both paths (--no-renames).
    left, _, _ = select(['tests/flow_discretization/petsc/old.cpp'])
    right, _, _ = select(['frontend/src/new.ts'])
    union, _, _ = select(['tests/flow_discretization/petsc/old.cpp', 'frontend/src/new.ts'])
    assert all(not (left[k] or right[k]) or union[k] for k in union)
    # Transitive local workflow limit, not merely direct calls (GitHub maximum 50).
    seen = set()
    def visit(path):
        for job in load(path).get('jobs', {}).values():
            used = job.get('uses', '')
            if used.startswith('./.github/workflows/') and used[2:] not in seen:
                seen.add(used[2:]); visit(used[2:])
    visit(router_path)
    assert len(seen) <= 50, ('too many unique reusable workflows', len(seen))
    print('WORKFLOW_MAP_OK', len(paths), 'entries; single automatic entry; reusable closure:', len(seen))
    print('ENTRY_INPUT_MATRIX_PERMISSIONS_AND_STEP_PARITY_OK')
if __name__ == '__main__':
    main()
'''
put('.github/ci/verify_workflows.py', verify)
put('.github/ci/README.md', '''# 单入口测试调度与旧入口映射

`workflow_map.json` 是逐 workflow 的审计清单：记录迁移前事件和路径/分支条件、完整手动输入及默认值、runner/平台矩阵、权限、environment、timeout、原 job 条件、中央 job 映射及语义指纹。不得为了匹配指纹而盲目更新清单；须先审计行为变化。

## 入口语义

- 常规 PR 只有 `pr_incremental_ci.yml` 自动运行。旧独立 PR job 进入这个 run，不再生成各自的空 workflow run。原先已由中央路由调用的 reusable workflow 不改执行体。
- 旧 `workflow_dispatch` 输入（包含 required/default/type/options）、手动签名/发布限制、平台矩阵、artifact、安装包验证、Clapeyron/ThermoPack 固定来源均保留。中央自动 job 不注入手动参数默认值，避免把手动签名/发布操作变成 PR 自动操作。
- 原来只有 PR 入口的门禁显式补参数为空的手动入口，不删除其唯一可执行路径。中央手动入口维持既有核心套件集合；特殊手动验证仍使用对应命名 workflow。
- 不把各旧 workflow 全部变成新增 reusable 调用，避免超过 GitHub 对单调用树唯一 reusable workflow 数量的限制；保留 legacy job 的 needs/output、defaults/env 和每 job 权限。

## 增量与失败

`plan.py` 仅采用同 PR、同 base、祖先可达且中央新版本完整成功的 checkpoint。失败、取消或排队的提交不前移基线；无证据、首次运行、Ready-for-review 回退累计 PR 差异。删除/重命名按旧新路径共同选测。未知可执行路径明确失败，不能静默漏测。首次迁移没有新版本 checkpoint，会验证累计受影响范围，不能为节省本轮时间伪造成功基线。

不受信任 fork 不在私人 runner 执行。Linux 科学测试保留 `mpmc_hnu`，Windows/macOS 保留官方平台。`Required CI result` 汇总失败/取消，不将 skipped 当作已跑过测试。分支保护不在此次修改范围内。

## 后续编写

先读 `.github/AGENTS.md` 与 `tests/AGENTS.md`；新增测试接现有 owner Gate，同步依赖规则、入口映射与反例。运行 `python3 .github/ci/verify_workflows.py`（PyYAML 6.0.2）。手动参数或矩阵变更必须更新映射并验证特殊路径。此次入口迁移不等同于全部测试计算已经去重；SW92 跨 workflow 的重复计算仍须单独审计，不能以单入口验收替代。
''')
subprocess.run(['python3', '.github/ci/verify_workflows.py'], check=True)
subprocess.run(['git', 'diff', '--check'], check=True)
print('MAPPED_WORKFLOWS', len(registry['workflows']))
print('LEGACY_CENTRAL_JOBS', len(legacy_jobs))
print('MANUAL_PARAMETERS', json.dumps({p: s['manual_inputs'] for p,s in registry['workflows'].items() if s['manual_inputs']}, ensure_ascii=False))
print('REPAIRED_PR_ONLY_ENTRIES', json.dumps([p for p,s in registry['workflows'].items() if 'entry_repair' in s]))
# Validate YAML semantics with the actual GitHub Actions linter if installed.
linter = Path.home() / 'go/bin/actionlint'
if linter.exists():
    Path('/tmp/actionlint.yaml').write_text('self-hosted-runner:\n  labels: [mpmc_hnu]\n')
    subprocess.run([str(linter), '-config-file', '/tmp/actionlint.yaml', '-shellcheck=', '-pyflakes='], check=True)
else:
    raise RuntimeError('actionlint is required before publishing this candidate')
changed = git('status', '--porcelain').splitlines()
assert all(line[3:].startswith('.github/') for line in changed), changed
# Create immutable candidate objects only; authorized connector advances the PR.
def api(path, payload):
    req = urllib.request.Request('https://api.github.com/repos/' + os.environ['GITHUB_REPOSITORY'] + path,
        data=json.dumps(payload).encode(), method='POST', headers={'Authorization': 'Bearer ' + os.environ['GH_TOKEN'], 'Accept': 'application/vnd.github+json', 'Content-Type': 'application/json'})
    with urllib.request.urlopen(req, timeout=30) as response:
        return json.load(response)
subprocess.run(['git', 'add', '.github'], check=True)
files = git('diff', '--cached', '--name-only').splitlines()
entries = [{'path': p, 'mode': '100644', 'type': 'blob', 'content': Path(p).read_text(encoding='utf-8')} for p in files]
tree = api('/git/trees', {'base_tree': git('rev-parse', BASE + '^{tree}'), 'tree': entries})
commit = api('/git/commits', {'message': 'refactor(ci): preserve legacy gate semantics in a single PR entry', 'tree': tree['sha'], 'parents': [BASE]})
print('CANDIDATE_COMMIT=' + commit['sha'])
print('CANDIDATE_TREE=' + tree['sha'])
print('CHANGED_PATHS=' + json.dumps(files))
