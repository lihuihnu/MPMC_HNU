"""One-shot, read-only-ref CI migration. Publish a candidate object, never a ref."""
import copy
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import urllib.request
import yaml

BASE = 'f75f703fe9203d3a9aef7e78a195c8da46f13366'
ROUTER = '.github/workflows/pr_incremental_ci.yml'
ROOT = Path('.')
def git(*args):
    return subprocess.check_output(['git', *args], text=True)
def load(text):
    value = yaml.safe_load(text)
    if isinstance(value, dict) and True in value:
        value['on'] = value.pop(True)
    return value
def dump(value):
    return yaml.safe_dump(value, sort_keys=False, allow_unicode=True, width=120)
def original(path):
    return git('show', BASE + ':' + path)
paths = git('ls-tree', '-r', '--name-only', BASE, '.github/workflows').splitlines()
originals = {p: original(p) for p in paths if p.endswith(('.yml', '.yaml'))}
workflows = {p: load(t) for p, t in originals.items()}
assert len(workflows) >= 60
router = copy.deepcopy(workflows[ROUTER])
old_jobs = router['jobs']
route_step = next(s for s in old_jobs['impact']['steps'] if s.get('id') == 'route')
rule_source = route_step['run'].split('event_name = os.environ.get("EVENT_NAME", "")')[0]
assert 'def route(' in rule_source and 'Thermodynamics selector regressions' in rule_source
rules = {}
exec(compile(rule_source, 'original-impact-rules', 'exec'), rules)
changes = {'.github/ci/impact_rules.py': rule_source}
registry = {'extra': {}, 'sources': {}, 'source_owner': {}, 'sw92': [], 'external': {}}

# Retain the existing selection rules. Add all former standalone PR selectors.
for path, wf in workflows.items():
    if path == ROUTER:
        continue
    events = wf.get('on') or {}
    assert isinstance(events, dict), (path, events)
    if 'pull_request' in events:
        event = events['pull_request'] or {}
        key = 'auto_' + Path(path).stem
        registry['extra'][key] = {'paths': event.get('paths', ['**']), 'ignore': event.get('paths-ignore', [])}
    for event in ('pull_request', 'pull_request_target', 'push', 'workflow_run', 'schedule'):
        events.pop(event, None)
    assert events, ('no manual/reusable entry remains', path)
    wf['on'] = events

# Source-level recipes are extracted without editing CMake/CTest commands.
def recipe(wf):
    jobs = [j for j in wf.get('jobs', {}).values() if 'steps' in j and 'strategy' in j]
    if not jobs:
        return None
    assert len(jobs) == 1, 'ambiguous SW92 matrix ownership'
    job = jobs[0]
    groups, extras, current = {}, [], None
    for step in job['steps']:
        command = step.get('run', '')
        matches = re.findall(r'cmake\s+-S\s+(tests/[\w/]+)\s+-B\s+([\w/.-]+)', command)
        if matches:
            assert len(matches) == 1, step
            current, build = matches[0]
            groups.setdefault(current, []).append(copy.deepcopy(step))
        elif re.search(r'cmake\s+--build|ctest\s+--test-dir', command):
            assert current, step
            groups[current].append(copy.deepcopy(step))
        elif 'actions/checkout@' in step.get('uses', '') or step.get('name') == 'Select Python interpreter':
            continue
        else:
            extras.append(copy.deepcopy(step))
    return job, groups, extras

# Fold the full routed Profile-C chain into one platform matrix, not just C2b1/C2b2.
sw_map = {}
for key, job in old_jobs.items():
    if key.startswith('sw92-') and 'uses' in job:
        path = job['uses'].removeprefix('./')
        if Path(path).name == '_sw92_phase_assigned_topology.yml':
            continue
        sw_map[key.replace('-', '_')] = path
for key in registry['extra']:
    name = key.removeprefix('auto_')
    if name.startswith('sw92_'):
        sw_map[key] = '.github/workflows/' + name + '.yml'
shared_path = '.github/workflows/_sw92_phase_assigned_topology.yml'
shared_job, shared_groups, shared_extras = recipe(workflows[shared_path])
for key, source in [
    ('sw92_phase_assigned_three_phase', 'tests/flash/sw92_phase_assigned_three_phase'),
    ('sw92_phase_assigned_three_phase_closure', 'tests/flash/sw92_phase_assigned_three_phase_closure'),
    ('sw92_phase_assigned_boundary', 'tests/flash/sw92_phase_assigned_boundary')]:
    registry['sources'][key] = {
        'sw92_phase_assigned_three_phase': ['tests/flash/sw92_phase_assigned_joint', 'tests/flash/sw92_phase_assigned_h_side_witness', source],
        'sw92_phase_assigned_three_phase_closure': ['tests/flash/sw92_phase_assigned_joint', 'tests/flash/sw92_phase_assigned_h_side_witness', 'tests/flash/sw92_phase_assigned_three_phase', source],
        'sw92_phase_assigned_boundary': ['tests/flash/sw92_phase_assigned_joint', 'tests/flash/sw92_phase_assigned_h_side_witness', 'tests/flash/sw92_phase_assigned_three_phase_closure', 'tests/flash/sw92_phase_assigned_no_w', 'tests/flash/sw92_phase_assigned_pt', source],
    }[key]
    registry['source_owner'][source] = key
    registry['sw92'].append(key)
leaf_steps = {s: shared_groups[s] for s in registry['source_owner']}
extra_steps = []
for key in registry['sw92']:
    extra_steps.extend((key, s) for s in shared_extras)
base_matrix = shared_job['strategy']['matrix']
for key, path in sw_map.items():
    parsed = recipe(workflows[path])
    assert parsed is not None, path
    job, groups, extras = parsed
    assert job['strategy']['matrix'] == base_matrix, ('SW92 platform matrix differs; preserve separately', path)
    assert groups, path
    source = next(iter(groups))
    assert source not in leaf_steps, ('duplicate primary ownership', path, source)
    leaf_steps[source] = groups[source]
    registry['source_owner'][source] = key
    registry['sources'][key] = list(groups)
    registry['sw92'].append(key)
    extra_steps.extend((key, s) for s in extras)
registry['external']['tests/flash/pt_split'] = 'pt_split'
for sources in registry['sources'].values():
    for source in sources:
        assert source in registry['source_owner'] or source in registry['external'], ('unowned dependency', source)

sw_job = copy.deepcopy(shared_job)
sw_job['needs'] = 'impact'
sw_job['if'] = "needs.impact.outputs.sw92_has_work == 'true'"
sw_job['name'] = 'SW92 affected gates / ${{ matrix.name }}'
sw_job['timeout-minutes'] = 90
sw_job['steps'] = [copy.deepcopy(s) for s in shared_job['steps'] if 'actions/checkout@' in s.get('uses', '') or s.get('name') == 'Select Python interpreter']
# Reference commands are unioned by exact command and execution environment.
unique_extras = {}
for origin, step in extra_steps:
    if 'inputs.run_' in str(step.get('if', '')):
        step.pop('if', None)
    command = step.get('run', '')
    lines = [line for line in command.splitlines() if line.strip()]
    pieces = lines if lines and all(re.match(r'\s*(?:"\$MPMC_PYTHON"|python3?)\s+tests/', line) for line in lines) else [None]
    for line in pieces:
        part = copy.deepcopy(step)
        if line is not None:
            part['run'] = line + '\n'
            part['name'] = 'Reference / ' + line.split()[1]
        signature = json.dumps({k: v for k, v in part.items() if k != 'name'}, sort_keys=True)
        if signature not in unique_extras:
            unique_extras[signature] = [part, set()]
        unique_extras[signature][1].add(origin)
for part, origins in unique_extras.values():
    condition = ' || '.join("contains(fromJSON(needs.impact.outputs.sw92_origins), '" + k + "')" for k in sorted(origins))
    old_if = part.pop('if', None)
    if old_if:
        condition = '(' + condition + ') && (' + str(old_if).removeprefix('${{').removesuffix('}}').strip() + ')'
    part['if'] = condition
    sw_job['steps'].append(part)
for source, steps in leaf_steps.items():
    assert sum('ctest ' in s.get('run', '') for s in steps) == 1, ('non-atomic test gate', source)
    for step in steps:
        step = copy.deepcopy(step)
        old_if = step.pop('if', '')
        condition = "contains(fromJSON(needs.impact.outputs.sw92_sources), '" + source + "')"
        if old_if and 'inputs.run_' not in str(old_if):
            condition += ' && (' + str(old_if).removeprefix('${{').removesuffix('}}').strip() + ')'
        step['if'] = condition
        step['name'] = source + ' / ' + step.get('name', 'validation')
        sw_job['steps'].append(step)

new_jobs = {k: copy.deepcopy(v) for k, v in old_jobs.items() if k != 'impact' and not k.startswith('sw92-')}
new_jobs['sw92'] = sw_job
# Inline former independent PR jobs. No new reusable-workflow nesting limit.
for path, wf in workflows.items():
    key = 'auto_' + Path(path).stem
    if key not in registry['extra'] or key in sw_map:
        continue
    original_wf = load(originals[path])
    jobs = copy.deepcopy(original_wf.get('jobs', {}))
    removed = {k for k, j in jobs.items() if j.get('uses', '').endswith('/_incremental_path_gate.yml')}
    rename = {k: key + '__' + k for k in jobs if k not in removed}
    dispatch = (original_wf.get('on', {}).get('workflow_dispatch') or {}).get('inputs', {})
    for old_key, job in jobs.items():
        if old_key in removed:
            continue
        dependencies = job.get('needs', [])
        if isinstance(dependencies, str):
            dependencies = [dependencies]
        job['needs'] = ['impact'] + [rename[n] for n in dependencies if n not in removed]
        serialized = json.dumps(job)
        for old, new in rename.items():
            serialized = re.sub(r'needs\.' + re.escape(old) + r'(?=[.\s\[])', 'needs.' + new, serialized)
        for old in removed:
            serialized = re.sub(r"needs\." + re.escape(old) + r"\.outputs\.affected\s*==\s*'true'", 'true', serialized)
        for name, spec in dispatch.items():
            literal = "fromJSON('" + json.dumps(spec.get('default', '')).replace("'", "''") + "')"
            serialized = serialized.replace('inputs.' + name, literal)
        job = json.loads(serialized)
        previous_if = str(job.pop('if', 'true')).removeprefix('${{').removesuffix('}}').strip()
        job['if'] = "needs.impact.outputs." + key + " == 'true' && (" + previous_if + ')'
        job['name'] = original_wf.get('name', key) + ' / ' + str(job.get('name', old_key))
        if 'uses' not in job:
            job['env'] = {**original_wf.get('env', {}), **job.get('env', {})}
            if 'defaults' in original_wf and 'defaults' not in job:
                job['defaults'] = copy.deepcopy(original_wf['defaults'])
        if 'permissions' in original_wf and 'permissions' not in job:
            job['permissions'] = copy.deepcopy(original_wf['permissions'])
        new_jobs[rename[old_key]] = job

changes['.github/ci/suites.json'] = json.dumps(registry, ensure_ascii=False, indent=2) + '\n'
plan_source = r'''"""Select only affected suites since the last successful compatible CI checkpoint."""
import json
import os
from pathlib import Path
import re
import runpy
import subprocess
import urllib.request

RULES = runpy.run_path('.github/ci/impact_rules.py')
REGISTRY = json.loads(Path('.github/ci/suites.json').read_text())
SCHEMA = 'single-entry-v1'
def affected(paths):
    result = RULES['route'](paths)
    for key, spec in REGISTRY['extra'].items():
        result[key] = any(any(RULES['compile_pattern'](p).match(path) for p in spec['paths']) and not any(RULES['compile_pattern'](p).match(path) for p in spec['ignore']) for path in paths)
    if any(path.startswith('.github/ci/') for path in paths):
        result = {key: True for key in result}
    elif paths and all(path in ('AGENTS.md', 'README.md') for path in paths):
        result = {key: False for key in result}
    origins, sources = set(), set()
    pending = [key for key in REGISTRY['sw92'] if result.get(key, False)]
    while pending:
        key = pending.pop()
        if key in origins:
            continue
        origins.add(key)
        for source in REGISTRY['sources'][key]:
            if source in REGISTRY['external']:
                result[REGISTRY['external'][source]] = True
            else:
                sources.add(source)
                pending.append(REGISTRY['source_owner'][source])
    unknown = [p for p in paths if p.startswith(('modules/', 'tests/', 'frontend/', 'api/')) and not p.endswith(('.md', '.txt'))]
    if unknown and not any(result.values()):
        raise RuntimeError('Unowned executable change; register its CI dependency closure: ' + repr(unknown))
    return result, sorted(origins), sorted(sources)
def git(*args):
    return subprocess.check_output(['git', *args], text=True).strip()
def ancestor(before, after):
    return subprocess.run(['git', 'merge-base', '--is-ancestor', before, after], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0
def api(path):
    request = urllib.request.Request('https://api.github.com/repos/' + os.environ['GITHUB_REPOSITORY'] + path, headers={'Authorization': 'Bearer ' + os.environ['GH_TOKEN'], 'Accept': 'application/vnd.github+json'})
    with urllib.request.urlopen(request, timeout=30) as response:
        return json.load(response)
def main():
    event = json.loads(Path(os.environ['GITHUB_EVENT_PATH']).read_text())
    pr = event.get('pull_request')
    head = pr['head']['sha'] if pr else git('rev-parse', 'HEAD')
    base = pr['base']['sha'] if pr else None
    before, mode = base, 'PR cumulative validation'
    if pr and event.get('action') == 'synchronize':
        try:
            runs = api('/actions/workflows/pr_incremental_ci.yml/runs?per_page=100&event=pull_request&status=success')['workflow_runs']
            for run in runs:
                if run['name'] != 'CI · affected tests' or run['head_branch'] != pr['head']['ref']:
                    continue
                if not any(p['number'] == pr['number'] and p['base']['sha'] == base for p in run.get('pull_requests', [])):
                    continue
                candidate = run['head_sha']
                if ancestor(base, candidate) and ancestor(candidate, head):
                    before, mode = candidate, 'delta since last successful CI checkpoint'
                    break
        except Exception as error:
            print('Checkpoint unavailable; use cumulative validation:', type(error).__name__)
    if before:
        raw = subprocess.check_output(['git', 'diff', '--name-only', '--no-renames', '-z', before, head])
        paths = [p for p in raw.decode().split('\0') if p]
        result, origins, sources = affected(paths)
        selected_ad = RULES['ad_suites_for'](paths)
        selected_thermo = RULES['thermo_suites_for'](paths) if result.get('thermodynamics_contracts') else []
        if any(p.startswith('.github/ci/') for p in paths):
            selected_ad = RULES['all_ad_suites']
            selected_thermo = RULES['all_thermo_suites']
    else:
        paths = ['.github/ci/manual-full-validation']
        result, origins, sources = affected(paths)
        selected_ad = RULES['all_ad_suites']
        selected_thermo = RULES['all_thermo_suites']
        mode = 'explicit manual full validation'
    if result.get('thermodynamics_contracts') and not selected_thermo:
        selected_thermo = RULES['all_thermo_suites']
    output = {key: 'true' if value else 'false' for key, value in result.items()}
    output.update(ad_has_work='true' if selected_ad else 'false', ad_matrix=json.dumps(RULES['ad_matrix'](selected_ad)), thermo_has_work='true' if selected_thermo else 'false', thermo_matrix=json.dumps(RULES['thermo_matrix'](selected_thermo)), sw92_has_work='true' if sources else 'false', sw92_origins=json.dumps(origins), sw92_sources=json.dumps(sources))
    with open(os.environ['GITHUB_OUTPUT'], 'a') as stream:
        for key, value in output.items():
            stream.write(key + '=' + value + '\n')
    selected = sorted(key for key, value in result.items() if value)
    report = {'schema': SCHEMA, 'mode': mode, 'base': before, 'head': head, 'selected': selected, 'sw92_sources_once': sources, 'changed_paths': paths}
    print(json.dumps(report, indent=2))
    with open(os.environ['GITHUB_STEP_SUMMARY'], 'a') as stream:
        stream.write('## CI impact and coverage\n```json\n' + json.dumps(report, indent=2) + '\n```\n')
if __name__ == '__main__':
    main()
'''
changes['.github/ci/plan.py'] = plan_source
outputs = copy.deepcopy(old_jobs['impact']['outputs'])
for key in list(registry['extra']) + ['sw92_has_work', 'sw92_origins', 'sw92_sources']:
    outputs[key] = '${{ steps.route.outputs.' + key + ' }}'
impact = {
    'name': 'Impact / governance', 'runs-on': 'ubuntu-24.04', 'timeout-minutes': 10,
    'if': "github.event_name != 'pull_request' || github.event.pull_request.head.repo.full_name == github.repository",
    'outputs': outputs,
    'steps': [
        {'uses': 'actions/checkout@3d3c42e5aac5ba805825da76410c181273ba90b1', 'with': {'persist-credentials': False, 'fetch-depth': 0}},
        {'name': 'Install pinned YAML validator', 'run': 'python3 -m pip install --disable-pip-version-check PyYAML==6.0.2'},
        {'name': 'Validate CI ownership and selection regressions', 'run': 'python3 .github/ci/test_governance.py'},
        {'name': 'Select affected dependency closure', 'id': 'route', 'env': {'GH_TOKEN': '${{ github.token }}'}, 'run': 'python3 .github/ci/plan.py'}]}
router['name'] = 'CI · affected tests'
router['on'] = {'pull_request': {'branches': ['main'], 'types': ['opened', 'synchronize', 'reopened', 'ready_for_review']}, 'workflow_dispatch': None}
router['permissions'] = {'contents': 'read', 'actions': 'read'}
router['jobs'] = {'impact': impact, **new_jobs}
router['jobs']['result'] = {
    'name': 'CI result', 'if': '${{ always() }}', 'needs': list(router['jobs']),
    'runs-on': 'ubuntu-24.04', 'timeout-minutes': 5,
    'steps': [{'name': 'Reject failed or cancelled validation', 'env': {'RESULTS': '${{ toJSON(needs) }}'}, 'shell': 'python', 'run': "import json, os\nr = json.loads(os.environ['RESULTS'])\nbad = {k: v['result'] for k, v in r.items() if v['result'] not in ('success', 'skipped')}\nif r['impact']['result'] != 'success':\n    bad['impact'] = r['impact']['result']\nprint(json.dumps(r, indent=2))\nif bad:\n    raise SystemExit('CI incomplete/failed: ' + repr(bad))\nprint('Selected gates completed; unselected gates were not executed.')\n"}]}
workflows[ROUTER] = router
for path, wf in workflows.items():
    text = dump(wf)
    if load(originals[path]) != wf:
        changes[path] = text

test_source = r'''"""CI invariants: a single trigger, owned dependency union, and valid DAGs."""
import importlib.util
import json
from pathlib import Path
import unittest
import yaml

def load(path):
    value = yaml.safe_load(path.read_text())
    if True in value:
        value['on'] = value.pop(True)
    return value
spec = importlib.util.spec_from_file_location('ci_plan', '.github/ci/plan.py')
plan = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plan)
class Governance(unittest.TestCase):
    def test_one_automatic_entry(self):
        entries = []
        for path in Path('.github/workflows').glob('*.yml'):
            events = load(path).get('on', {})
            if any(e in events for e in ('pull_request', 'pull_request_target', 'push', 'schedule', 'workflow_run')):
                entries.append(path.name)
        self.assertEqual(entries, ['pr_incremental_ci.yml'])
    def test_dag_and_result(self):
        jobs = load(Path('.github/workflows/pr_incremental_ci.yml'))['jobs']
        for key, job in jobs.items():
            needs = job.get('needs', [])
            if isinstance(needs, str): needs = [needs]
            self.assertTrue(set(needs) <= jobs.keys(), key)
            self.assertNotIn(key, needs)
        self.assertEqual(set(jobs['result']['needs']), set(jobs) - {'result'})
    def test_sw92_exactly_once(self):
        job = load(Path('.github/workflows/pr_incremental_ci.yml'))['jobs']['sw92']
        commands = [s['run'] for s in job['steps'] if 'ctest ' in s.get('run', '')]
        self.assertEqual(len(commands), len(set(commands)))
        self.assertEqual(len(commands), len(plan.REGISTRY['source_owner']))
        for key in plan.REGISTRY['sw92']:
            paths = ['.github/ci/closure-regression']
            result, origins, sources = plan.affected(paths)
            self.assertEqual(len(sources), len(set(sources)))
            self.assertEqual(set(sources), set(plan.REGISTRY['source_owner']))
    def test_document_only(self):
        result, origins, sources = plan.affected(['AGENTS.md', 'README.md'])
        self.assertFalse(any(result.values()))
        self.assertFalse(origins or sources)
    def test_shared_sw92_dependency(self):
        result, origins, sources = plan.affected(['tests/support/sw92/test_support.hpp'])
        self.assertIn('tests/flash/sw92_phase_assigned_joint', sources)
        self.assertIn('tests/flash/sw92_phase_assigned_pt', sources)
        self.assertIn('tests/flash/sw92_profile_c_phase_set', sources)
        self.assertEqual(len(sources), len(set(sources)))
    def test_platforms_preserved(self):
        matrix = load(Path('.github/workflows/pr_incremental_ci.yml'))['jobs']['sw92']['strategy']['matrix']['include']
        self.assertEqual([m['runner'] for m in matrix], ['mpmc_hnu', 'mpmc_hnu', 'windows-2022'])
        self.assertEqual([m['sanitizer'] for m in matrix], ['ON', 'OFF', 'OFF'])
if __name__ == '__main__':
    unittest.main(verbosity=2)
'''
changes['.github/ci/test_governance.py'] = test_source
policy = '''

## 12. 测试编写与 CI 治理强制规范

本节取代前文中允许多个独立 PR workflow 自动触发的旧约定。

- **单一自动入口。** 只有 `.github/workflows/pr_incremental_ci.yml` 监听 PR；新增测试默认接入已有 CTest target/标签与中央路由，不再新建独立 `pull_request/push/workflow_run/schedule` 入口。已有独立 workflow 只保留明确的手动或 reusable 入口，不创建先启动再 skipped 的无关 workflow runs。
- **按影响和证据选测。** 正常增量以同一 PR、同一 base 且可达的最近成功 CI checkpoint 为起点；失败、取消和仍运行的提交不能充当已验证基线。首次验证、无法确认 checkpoint、重开和 Ready-for-review 使用累计 PR diff。删除、重命名、公共头、CMake、共享 fixture、oracle、协议和适配层的必要下游必须纳入。未知可执行改动不得静默跳过。
- **唯一执行 ownership。** 每个 `(测试 Gate, 编译器/平台, 配置, sanitizer/其他验证模式)` 在一次自动 CI 内只有一个 owner。依赖通过 Gate 集合并集表达，不允许上层 workflow 重复 configure/build/ctest 下层 Gate。不同验证模式不得为了去重合并；软件结构测试不代替独立物理验证。共享 oracle 同一配置按相同命令去重，来源与参数不变。
- **新增测试必须同步路由回归。** 新增或修改 CTest 时明确测试不变量、源目录、目标、测试名/标签、上游输入与必要下游；同步 `.github/ci/impact_rules.py` / `suites.json`，并提供直接命中、共享依赖、无关提交、删除/重命名和重复 closure 的选择回归。公共头自包含测试应是独立 CTest，不使用只编译却未执行的伪覆盖。测试依赖不得传播给库使用者。
- **最小充分而非少测。** 保留科学容差、golden/reference、编译告警、GCC sanitizer、Clang 与 MSVC 的适用验证。不得通过删测试、降低精度、改 golden 或忽略失败来减少成本。文档/治理规则变更只执行相关结构检查；修改执行命令或公共依赖时验证受影响矩阵。
- **runner ownership。** Linux 编译和数值测试使用私有 `mpmc_hnu`；Windows/macOS 使用官方 runner。官方 Linux 只承担不编译科学代码的路径分析、结构回归和结果汇总。来自不受信任 fork 的代码不得自动进入私有 runner。缓存必须包含工具链、依赖、配置及 sanitizer 的身份，不复用不相容构建。
- **结果语义。** Actions success 不自动等同于 CTest 已执行。汇报 selected、unselected、实际完成、失败、取消和未验证状态，并记录 diff 起止 SHA。汇总门禁不能把失败/取消当作通过；必需检查只绑定稳定的 `CI result`，仓库保护设置由负责人单独授权管理，不在本次治理中修改。
- **一批提交、一轮验收。** 写入前完成静态审计与路由自测，多文件统一提交，不用逐文件 push 反复取消 runner。预检候选与正式科学测试是不同证据，预检不能冒充正式门禁通过。

CI 自检入口：`python3 -m pip install PyYAML==6.0.2`，然后 `python3 .github/ci/test_governance.py`。手动完整验证使用中央 workflow 的 `workflow_dispatch`；保留的专项手动入口不用于常规提交重复验证。
'''
changes['AGENTS.md'] = original('AGENTS.md') + policy
assert all(p.startswith('.github/') or p == 'AGENTS.md' for p in changes)
for path, text in changes.items():
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text, encoding='utf-8')
subprocess.run(['python3', '.github/ci/test_governance.py'], check=True)
# Recheck source inventory and original scientific commands before publication.
assert all('ctest ' in ''.join(s.get('run', '') for s in steps) for steps in leaf_steps.values())
print('PREVIEW:', json.dumps({'changed_files': sorted(changes), 'automatic_workflows': 1, 'sw92_unique_sources': sorted(leaf_steps), 'independent_selectors': len(registry['extra'])}, indent=2))
repo = os.environ['GITHUB_REPOSITORY']
event = json.loads(Path(os.environ['GITHUB_EVENT_PATH']).read_text())
parent = event['pull_request']['head']['sha']
def post(endpoint, payload):
    request = urllib.request.Request('https://api.github.com/repos/' + repo + endpoint, data=json.dumps(payload).encode(), headers={'Authorization': 'Bearer ' + os.environ['GH_TOKEN'], 'Accept': 'application/vnd.github+json', 'Content-Type': 'application/json'}, method='POST')
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)
elements = [{'path': p, 'mode': '100644', 'type': 'blob', 'content': text} for p, text in changes.items()]
tree = post('/git/trees', {'base_tree': git('rev-parse', BASE + '^{tree}'), 'tree': elements})
commit = post('/git/commits', {'message': 'refactor(ci): enforce single-entry affected-test ownership\n\nPreserve scientific commands and validation modes; dedupe SW92 closures.\nAdd checkpoint-aware selection, stable result gate and mandatory test rules.\nCandidate metadata tests passed; full platform validation is pending.', 'tree': tree['sha'], 'parents': [parent]})
print('CI_GOVERNANCE_CANDIDATE=' + commit['sha'])
print('Candidate only: no branch ref was changed and no scientific test pass is claimed.')
