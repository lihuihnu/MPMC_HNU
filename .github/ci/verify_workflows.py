"""Executable entry/mapping guard. Update the map only with an audited CI change."""
import hashlib
import json
from pathlib import Path
import runpy
import subprocess
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
            if spec.get('retained_blob_sha'):
                actual_blob = subprocess.check_output(['git', 'hash-object', path], text=True).strip()
                assert actual_blob == spec['retained_blob_sha'], ('workflow blob differs from audited map', path)
            else:
                assert digest(wf) == spec['retained_hash'], ('workflow semantics differ from audited map', path)
            if 'workflow_dispatch' in spec['original_events']:
                assert events['workflow_dispatch'] == spec['original_events']['workflow_dispatch'], path
        for job in spec.get('central_hashes', {}):
            assert job in root['jobs'], ('mapped central job missing', job)
    assert auto == [router_path], ('multiple automatic workflow entries', auto)
    # Private Linux runners are reserved for audited long-running/stable-hardware gates.
    private_allowlist = {
        ('.github/workflows/cpa_performance_audit.yml', 'paired-baseline'),
        (router_path, 'legacy_cpa_performance_audit__paired-baseline'),
    }
    private_seen = set()
    for workflow_path in sorted(paths):
        wf = load(workflow_path)
        for job_name, job in (wf.get('jobs') or {}).items():
            serialized = json.dumps(job, sort_keys=True)
            uses_private = 'mpmc_hnu' in serialized or 'self-hosted' in serialized
            if uses_private:
                pair = (workflow_path, job_name)
                assert pair in private_allowlist, ('private runner used outside long-test allowlist', pair)
                private_seen.add(pair)
    assert private_seen == private_allowlist, ('private runner allowlist drift', private_seen, private_allowlist)

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
    # Central workflow-level cancellation supersedes old heads; inline matrix jobs
    # must not carry copied concurrency groups that cancel siblings in the same run.
    for name, job in root['jobs'].items():
        if name.startswith('legacy_') and 'strategy' in job:
            assert 'concurrency' not in job, ('matrix job may self-cancel siblings', name)

    # Direct consumers of selected_phase_density must both be selected.
    impact_rules = runpy.run_path('.github/ci/impact_rules.py')
    density = impact_rules['route']([
        'modules/thermodynamics/include/mpmc/thermodynamics/selected_phase_density.hpp'
    ])
    assert density['flow_core'] and density['flow_discretization_petsc']

    # Profile-C PT is owned by the shared topology closure, not a second reusable job.
    assert 'sw92-profile-c-pt' not in root['jobs']
    topology = root['jobs']['sw92-phase-assigned-three-phase-topology']
    assert 'sw92_profile_c_pt' in str(topology.get('if', ''))
    assert 'sw92_profile_c_pt' in json.dumps(topology.get('with', {}))

    # Result gate must compare selected outputs, not merely accept skipped jobs.
    result_text = json.dumps(root['jobs']['result'], ensure_ascii=False)
    assert 'IMPACT_JSON' in result_text
    assert 'verify_result.py' in result_text
    result_guard = runpy.run_path('.github/ci/verify_result.py')
    fake_catalog = {'workflows': {}}
    try:
        result_guard['validate'](
            {'impact': {'result': 'success'}, 'flow-core': {'result': 'skipped'}},
            {'trusted': 'true', 'flow_core': 'true'},
            fake_catalog,
        )
    except AssertionError:
        pass
    else:
        raise AssertionError('selected-but-skipped result was accepted')
    result_guard['validate'](
        {'impact': {'result': 'success'}, 'flow-core': {'result': 'success'}},
        {'trusted': 'true', 'flow_core': 'true'},
        fake_catalog,
    )

    # Routed reusable workflows publish explicit validated outputs. Caller-level
    # cancellation must not hide a successful validated matrix, while missing
    # validated evidence must still fail.
    result_guard['validate'](
        {
            'impact': {'result': 'success'},
            'ad': {'result': 'cancelled', 'outputs': {'validated': 'true'}},
            'pt-stability': {'result': 'cancelled', 'outputs': {'validated': 'true'}},
        },
        {'trusted': 'true', 'ad_has_work': 'true', 'pt_stability': 'true'},
        fake_catalog,
    )
    try:
        result_guard['validate'](
            {
                'impact': {'result': 'success'},
                'ad': {'result': 'cancelled', 'outputs': {}},
            },
            {'trusted': 'true', 'ad_has_work': 'true'},
            fake_catalog,
        )
    except AssertionError:
        pass
    else:
        raise AssertionError('routed reusable cancellation without validated output was accepted')

    # Central Profile-C calls skip only dependencies already owned by the topology closure.
    assert root['jobs']['sw92-profile-c-phase-set']['with']['dependencies_prevalidated'] is True
    assert root['jobs']['sw92-profile-c-sensitivity']['with']['dependencies_prevalidated'] is True
    assert root['jobs']['sw92-phase-assigned-no-w']['with']['dependencies_prevalidated'] is True

    print('WORKFLOW_MAP_OK', len(paths), 'entries; single automatic entry; reusable closure:', len(seen))
    print('ENTRY_INPUT_MATRIX_PERMISSIONS_AND_STEP_PARITY_OK')
if __name__ == '__main__':
    main()
