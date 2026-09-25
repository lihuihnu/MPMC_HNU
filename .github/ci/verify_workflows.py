"""Executable entry/mapping guard. Update the map only with an audited CI change."""
import hashlib
import json
from pathlib import Path
import re
import runpy
import subprocess
import tempfile
import yaml

def load(path):
    obj = yaml.safe_load(Path(path).read_text(encoding='utf-8'))
    if True in obj:
        obj['on'] = obj.pop(True)
    return obj

def digest(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False).encode()).hexdigest()

def job_run_text(job):
    runs = []
    for step in job.get('steps') or []:
        run = step.get('run')
        if run:
            runs.append(str(run))
    return '\n'.join(runs)

def workflow_run_text(path):
    workflow = load(path)
    return '\n'.join(
        job_run_text(job)
        for job in (workflow.get('jobs') or {}).values()
    )

def configured_ctest_names(source_dir):
    with tempfile.TemporaryDirectory(prefix='mpmc-ci-inventory-') as build_dir:
        configured = subprocess.run(
            ['cmake', '-S', source_dir, '-B', build_dir],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        assert configured.returncode == 0, (
            'CTest inventory configure failed',
            source_dir,
            configured.stdout,
        )
        listed = subprocess.check_output(
            ['ctest', '--test-dir', build_dir, '-N'],
            text=True,
        )
    names = re.findall(r'Test\s+#\d+:\s+(\S+)', listed)
    assert names, ('no registered CTests discovered', source_dir, listed)
    return sorted(set(names))

def assert_run_text_covers_registered_ctests(registered, run_text, owner, required_targets):
    for target in required_targets:
        assert target in run_text, (
            'registered required test target omitted from workflow build',
            owner,
            target,
        )

    raw_patterns = re.findall(r"-R\s+['\"]([^'\"]+)['\"]", run_text)
    assert raw_patterns, ('workflow has no CTest selection regex', owner)
    patterns = [re.compile(pattern) for pattern in raw_patterns]
    omitted = [
        name
        for name in registered
        if not any(pattern.search(name) for pattern in patterns)
    ]
    assert not omitted, (
        'registered CTests omitted by authoritative workflow selection',
        owner,
        omitted,
    )

def assert_workflow_covers_registered_ctests(source_dir, workflow_path, required_targets):
    registered = configured_ctest_names(source_dir)
    assert_run_text_covers_registered_ctests(
        registered,
        workflow_run_text(workflow_path),
        workflow_path,
        required_targets,
    )
    return registered

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
    # Private Linux runners are reserved for audited long-running gates and
    # explicitly authorized PETSc/MPI integration/solver gates.
    private_allowlist = {
        ('.github/workflows/cpa_performance_audit.yml', 'paired-baseline'),
        ('.github/workflows/flow_discretization_petsc.yml', 'distributed'),
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

    # Registered-test inventory must be closed by the authoritative build target
    # list and CTest regex. This catches the silent failure mode where CMake
    # registers a required test but a narrower workflow selection never builds
    # or runs it.
    discretization_core_tests = assert_workflow_covers_registered_ctests(
        'tests/discretization/core',
        '.github/workflows/discretization_core.yml',
        {'mpmc_discretization_core_all'},
    )
    expected_well_control_tests = {
        'well.core.single_well_control_policy.rate_hold_and_switch',
        'well.core.single_well_control_policy.rate_pressure_equality_boundary',
        'well.core.single_well_control_policy.capacity_deadband_and_equality',
        'well.core.single_well_control_policy.pressure_deadband_and_equality',
        'well.core.single_well_control_policy.disabled_reactivation',
        'well.core.single_well_control_policy.accepted_state',
        'well.core.single_well_control_policy.invalid_policy',
        'well.core.single_well_control_policy.header_self_contained',
    }
    assert expected_well_control_tests <= set(discretization_core_tests), (
        'single-well control-policy CTest registration drift',
        sorted(expected_well_control_tests - set(discretization_core_tests)),
    )
    assert_run_text_covers_registered_ctests(
        discretization_core_tests,
        job_run_text(root['jobs']['legacy_discretization_core__core']),
        'pr_incremental_ci.yml:legacy_discretization_core__core',
        {'mpmc_discretization_core_all'},
    )

    flow_core_tests = assert_workflow_covers_registered_ctests(
        'tests/flow/core',
        '.github/workflows/flow_core.yml',
        {
            'mpmc_flow_absent_phase_thermodynamics_tests',
            'mpmc_flow_pr76_selected_phase_property_tests',
            'mpmc_flow_sw92_selected_phase_property_tests',
        },
    )
    flow_discretization_tests = assert_workflow_covers_registered_ctests(
        'tests/flow_discretization/core',
        '.github/workflows/flow_discretization.yml',
        {'mpmc_flow_discretization_cell_source_tests'},
    )
    flow_core_run_text = workflow_run_text('.github/workflows/flow_core.yml')
    for token in (
        'tests/flow/core/reference_li_firoozabadi_sour_gas_flow.py',
        '--precision 80',
        '--precision 96',
    ):
        assert token in flow_core_run_text, (
            'Li-Firoozabadi independent flow oracle lost authoritative CI ownership',
            token,
        )

    # The flow-discretization PETSc suite is one monolithic executable.  Source
    # membership alone is not execution proof, so keep an explicit entry-point
    # guard for the SW92 transactional restart regression and public-header probe.
    petsc_cmake_text = Path('tests/flow_discretization/petsc/CMakeLists.txt').read_text(encoding='utf-8')
    petsc_entry_text = Path('tests/flow_discretization/petsc/distributed_component_conservation_test.cpp').read_text(encoding='utf-8')
    for token in (
        'sw92_transactional_phase_transition_restart_test.cpp',
        'sw92_transactional_phase_transition_restart_header.cpp',
    ):
        assert token in petsc_cmake_text, (
            'SW92 transactional restart source lost PETSc target ownership',
            token,
        )
    assert re.search(
        r'int\s+main\s*\([^)]*\)\s*\{[\s\S]*?sw92_transactional_phase_transition_restart_test\s*\(\s*\)\s*;',
        petsc_entry_text,
    ), 'SW92 transactional restart regression is compiled but not executed'
    assert re.search(
        r'void\s+headers\s*\(\s*\)\s*\{[\s\S]*?sw92_transactional_phase_transition_restart_header\s*\(\s*\)',
        petsc_entry_text,
    ), 'SW92 transactional restart header probe is compiled but not executed'

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
    results, _, _ = select(['tests/flow/core/reference_li_firoozabadi_sour_gas_flow.py'])
    assert results['flow_core'] and not results['flow_discretization']
    results, _, _ = select(['.github/workflows/flow_core.yml'])
    assert results['flow_core']
    results, _, _ = select(['.github/workflows/flow_discretization.yml'])
    assert results['flow_discretization']
    results, _, _ = select(['modules/flow/discretization/include/mpmc/flow_discretization/cell_source.hpp'])
    assert results['flow_discretization'] and not results['flow_core'] and not results['flow_discretization_petsc']
    results, _, _ = select(['modules/flow/discretization/petsc/include/mpmc/flow_discretization_petsc/physical_timestep_driver.hpp'])
    assert results['flow_discretization_petsc'] and not results['flow_core'] and not results['flow_discretization']
    results, _, _ = select(['modules/flow/discretization/petsc/include/mpmc/flow_discretization_petsc/cell_scoped_mixed_cardinality_evaluator_dispatcher.hpp'])
    assert results['flow_discretization_petsc'] and not results['flow_core'] and not results['flow_discretization']
    results, _, _ = select(['modules/mesh/petsc/include/mpmc/mesh_petsc/adapter.hpp'])
    assert results['legacy_mesh_petsc'] and not results['legacy_mesh_core']
    results, _, _ = select(['modules/discretization/petsc/include/mpmc/discretization_petsc/adapter.hpp'])
    assert results['legacy_mesh_petsc'] and results['flow_discretization_petsc'] and not results['legacy_discretization_core']
    results, _, _ = select(['modules/well/include/mpmc/well/peaceman_well_index_3d.hpp'])
    assert results['legacy_discretization_core'] and not results['legacy_mesh_petsc']
    results, _, _ = select(['tests/well/core/single_well_control_policy_test.cpp'])
    assert results['legacy_discretization_core'] and not results['flow_discretization']
    results, _, _ = select(['modules/well/discretization/include/mpmc/well_discretization/hydraulic_conductance.hpp'])
    assert results['flow_discretization'] and not results['legacy_discretization_core'] and not results['flow_core']
    results, _, _ = select(['tests/well/discretization/hydraulic_conductance_test.cpp'])
    assert results['flow_discretization'] and not results['legacy_discretization_core']
    results, _, _ = select(['modules/well/discretization/petsc/include/mpmc/well_discretization_petsc/fixed_bhp_well_source_evaluator.hpp'])
    assert results['flow_discretization_petsc'] and not results['flow_discretization'] and not results['flow_core']
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
    print(
        'REGISTERED_CTEST_INVENTORY_OK',
        'discretization-core=', len(discretization_core_tests),
        'flow-core=', len(flow_core_tests),
        'flow-discretization=', len(flow_discretization_tests),
    )
    print('ENTRY_INPUT_MATRIX_PERMISSIONS_AND_STEP_PARITY_OK')
if __name__ == '__main__':
    main()
