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

def assert_mesh_independent_readers(root, select):
    owner = 'legacy_mesh_external_compatibility'
    directory = Path('tests/mesh/external_compatibility')
    manual = load('.github/workflows/mesh_external_compatibility.yml')['jobs']['compatibility']
    central = root['jobs'][owner + '__compatibility']
    # Check execution, not just a path match or a compiled-but-unused producer.
    names = ('Install pinned independent mesh readers',
             'Independently read 2D exports with Gmsh and VTK',
             'Preserve independent reader evidence',
             'Independently read 3D exports with Gmsh and VTK',
             'Preserve 3D independent reader evidence',
             'Validate externally generated 3D mixed mesh chains',
             'Preserve 3D chain evidence',
             'Validate externally generated 2D mixed mesh chains',
             'Preserve 2D chain evidence',
             'Install pinned independent GRDECL reader',
             'Validate independent GRDECL complete chains',
             'Validate computational preparation through three formats',
             'Validate independent VTU named group chains',
             'Preserve GRDECL and computational evidence',
             'Prepare pinned optional HDF5 dependency',
             'Validate MRST face mesh HDF5 and polyhedral VTU',
             'Preserve MRST face mesh evidence')
    for name in names:
        left = [step for step in manual['steps'] if step.get('name') == name]
        right = [step for step in central['steps'] if step.get('name') == name]
        assert len(left) == len(right) == 1 and left == right, ('mesh reader step parity', name)
    entry = (directory / 'external_mesh_compatibility.cpp').read_text(encoding='utf-8')
    for dimension, marker in ((2, 'files=12 reports=6 negative_controls=5'),
                              (3, 'files=24 reports=12 negative_controls=10')):
        run = next(step for step in central['steps'] if step.get('name') ==
                   f'Independently read {dimension}D exports with Gmsh and VTK')
        assert not run.get('continue-on-error') and not central.get('continue-on-error')
        for token in ('set -euo pipefail', f'verify_{dimension}d_readers.py --producer',
                      '/mpmc_mesh_external_compatibility', '--output-dir', marker):
            assert token in run['run'], ('mesh reader failure/execution guard', token)
        assert f'emit_{dimension}d_exports.cpp' in (directory / 'CMakeLists.txt').read_text(encoding='utf-8')
        assert f'--emit-{dimension}d' in entry and f'emit_{dimension}d_exports(argv[2])' in entry
    for dimension, marker in ((2, 'inputs=2 chains=4 reports=4 disconnected_chains=4'),
                              (3, 'inputs=4 chains=8 reports=8 disconnected_chains=8')):
        chain = next(step for step in central['steps'] if step.get('name') ==
                     f'Validate externally generated {dimension}D mixed mesh chains')
        assert not chain.get('continue-on-error')
        for token in ('set -euo pipefail', f'verify_{dimension}d_chains.py --producer',
                      '/mpmc_mesh_external_compatibility', '--output-dir', marker):
            assert token in chain['run'], ('chain execution/failure guard', token)
        assert f'convert_{dimension}d_file.cpp' in (directory / 'CMakeLists.txt').read_text(encoding='utf-8')
        assert f'--convert-{dimension}d' in entry and f'convert_{dimension}d_file(argv[2], argv[3], argv[4])' in entry
    for filename in ('emit_2d_exports.cpp', 'verify_2d_readers.py', 'requirements-readers.txt',
                     'emit_3d_exports.cpp', 'verify_3d_readers.py',
                     'convert_3d_file.cpp', 'verify_3d_chains.py',
                     'convert_2d_file.cpp', 'verify_2d_chains.py',
                     'convert_grdecl_file.cpp', 'verify_grdecl_chains.py', 'requirements-grdecl-reader.txt',
                     'validate_computational_scale.cpp', 'convert_vtu_groups.cpp',
                     'verify_vtu_group_chains.py', 'vtu_groups.py',
                     'verify_face_mesh_bridge.py', 'mrst_folder_audit.py'):
        path = (directory / filename).as_posix()
        for action in ('opened', 'synchronize', 'ready_for_review'):
            chosen, _, _ = select([path], action=action)
            assert chosen[owner] and not chosen['legacy_mesh_core'], (path, action)
        before, _, _ = select([path])  # Deletion retains the old path.
        renamed, _, _ = select([path, (directory / ('renamed_' + filename)).as_posix()])
        assert all(not value or renamed[key] for key, value in before.items())
    for header in ('gmsh_4_1.hpp', 'vtu.hpp', 'linear_cell_mesh_2d.hpp', 'mesh_exchange_io.hpp',
                   'gmsh_4_1_3d.hpp', 'vtu_3d.hpp', 'linear_cell_mesh_3d.hpp', 'mesh_exchange.hpp'):
        assert select(['modules/mesh/include/mpmc/mesh/' + header])[0][owner], header
    assert not select(['modules/mesh/README.md'])[0][owner]
    requirements = (directory / 'requirements-readers.txt').read_text(encoding='utf-8')
    for line in requirements.splitlines():
        if line and not line.startswith('#'):
            assert re.fullmatch(r'[a-zA-Z0-9_-]+==[0-9][a-zA-Z0-9.]*', line), line
    assert 'gmsh==' in requirements and 'vtk==' in requirements
    assert 'h5py==3.15.1' in requirements and 'psutil==7.1.0' in requirements
    for job in (manual, central, root['jobs']['legacy_mesh_core__core']):
        prep=next(step for step in job['steps'] if step.get('name')=='Prepare pinned optional HDF5 dependency')
        assert 'tests/mesh/core/prepare_hdf5.py --root' in prep['run'] and not prep.get('continue-on-error')
        commands='\n'.join(step.get('run','') for step in job['steps'])
        assert '-DMPMC_MESH_WITH_HDF5=ON' in commands and '-DHDF5_USE_STATIC_LIBRARIES=ON' in commands
    bridge=next(step for step in central['steps'] if step.get('name')=='Validate MRST face mesh HDF5 and polyhedral VTU')
    for token in ('set -euo pipefail','verify_face_mesh_bridge.py --converter','negative_controls=16','grep -Fqx'):
        assert token in bridge['run'], ('face bridge execution guard',token)
    assert not bridge.get('continue-on-error')
    for path in ('modules/mesh/src/face_mesh_hdf5.cpp','modules/mesh/tools/mesh_convert.cpp',
                 'modules/mesh/matlab/export_mpmc_mesh.m','modules/mesh/include/mpmc/mesh/face_mesh.hpp',
                 'tests/mesh/core/prepare_hdf5.py'):
        gates=select([path])[0]
        assert gates[owner] and gates['legacy_mesh_core'],path
    for name, script, marker in (
            ('Validate independent GRDECL complete chains', 'verify_grdecl_chains.py',
             'chains=2 reports=2 negative_controls=13'),
            ('Validate independent VTU named group chains', 'verify_vtu_group_chains.py',
             'chains=3 reports=3 negative_controls=15')):
        step = next(step for step in central['steps'] if step.get('name') == name)
        assert not step.get('continue-on-error')
        for token in ('set -euo pipefail', script + ' --producer', '--output-dir', marker):
            assert token in step['run'], ('new independent chain execution/failure guard', token)
    scale = next(step for step in central['steps'] if step.get('name') ==
                 'Validate computational preparation through three formats')
    assert not scale.get('continue-on-error')
    for token in ('set -euo pipefail', 'for format in gmsh vtu grdecl', '--computational-scale', 'cells=64'):
        assert token in scale['run'], ('computational preparation scale guard', token)
    # Protect the transitive header closure of the only external producer.
    pending = list(directory.glob('*.cpp'))
    visited = set()
    while pending:
        path = pending.pop()
        if path in visited:
            continue
        visited.add(path)
        for header in re.findall(r'#include\s+[<"](mpmc/mesh/[^>"]+)[>"]', path.read_text(encoding='utf-8')):
            dependency = Path('modules/mesh/include') / header
            assert select([dependency.as_posix()])[0][owner], dependency
            pending.append(dependency)
    for line in (directory / 'requirements-grdecl-reader.txt').read_text(encoding='utf-8').splitlines():
        if line and not line.startswith('#'):
            assert re.fullmatch(r'[a-zA-Z0-9_-]+==[0-9][a-zA-Z0-9.]*', line), line

def assert_root_build_selection(select):
    direct = [
        'tests/build/root_libraries/CMakeLists.txt',
        'tests/build/root_libraries/consumer.cpp',
        'tests/build/root_libraries/spatial_consumer.cpp',
        'tests/build/root_libraries/flow_consumer.cpp',
        'tests/build/root_libraries/flow_discretization_consumer.cpp',
        'tests/build/root_libraries/well_consumer.cpp',
        'tests/build/root_libraries/well_discretization_consumer.cpp',
        'tests/build/root_libraries/verify.py',
        'modules/thermodynamics/CMakeLists.txt',
        'modules/flash/CMakeLists.txt',
        'modules/mesh/CMakeLists.txt',
        'modules/discretization/CMakeLists.txt',
        'modules/flow/CMakeLists.txt',
        'modules/flow/discretization/CMakeLists.txt',
        'modules/well/CMakeLists.txt',
        'modules/well/discretization/CMakeLists.txt',
        'tests/build/root_libraries/README.md',
    ]
    # Discover the probe's actual project-header closure, so a future include
    # cannot silently escape the explicit build-consumer ownership list.
    pending = [Path(path) for path in direct if path.endswith('.cpp')]
    visited = set()
    while pending:
        path = pending.pop()
        if path in visited:
            continue
        visited.add(path)
        text = path.read_text(encoding='utf-8')
        for header in re.findall(r'#include\s+[<"](mpmc/[^>"]+)[>"]', text):
            module = header.split('/')[1]
            # Bridge namespaces live under their owning module, not at root.
            module_dir = {
                'flow_discretization': 'flow/discretization',
                'well_discretization': 'well/discretization',
            }.get(module, module)
            pending.append(Path('modules') / module_dir / 'include' / header)
    direct.extend(str(path).replace('\\', '/') for path in visited)
    for action in ('opened', 'reopened', 'synchronize', 'ready_for_review'):
        for path in direct:
            assert 'arithmetic' in select([path], action=action)[1], (action, path)
        for path in ('modules/mesh/include/mpmc/mesh/vtu.hpp',
                     'modules/well/include/mpmc/well/single_well_control_policy.hpp',
                     'modules/well/discretization/petsc/CMakeLists.txt',
                     'modules/flow/include/mpmc/flow/sw92_co2_water_properties.hpp',
                     'modules/flash/include/mpmc/flash/sw92_profile_c_phase_set.hpp',
                     'tests/unknown/new.cpp', '.github/AGENTS.md'):
            assert not select([path], action=action)[1], ('unrelated AD fanout', path)
        # Includes old/new paths for deletion and rename selection; preserve the
        # science/legacy owners as well as the newly selected arithmetic suite.
        paths = ['tests/build/root_libraries/old.cpp', 'frontend/src/new.ts']
        left, right = (select([path], action=action) for path in paths)
        combined = select(paths, action=action)
        assert combined[1] == sorted(set(left[1]) | set(right[1]))
        assert all(combined[0][key] == (left[0][key] or right[0][key])
                   for key in combined[0])
    ad = load('.github/workflows/ad.yml')
    consumers = [step for job in ad['jobs'].values() for step in job.get('steps', [])
                 if 'tests/build/root_libraries/verify.py' in str(step.get('run', ''))]
    assert len(consumers) == 1 and consumers[0]['if'] == "matrix.suite == 'arithmetic'"

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

def assert_root_readme_selection(select):
    # Root overview edits have no scientific execution dependency. Topic docs,
    # source, fixtures and build inputs keep their existing owners even when
    # README.md is also changed (including either side of a rename).
    science_paths = (
        'modules/flash/include/mpmc/flash/sw92_profile_c_phase_set.hpp',
        'tests/flash/sw92_profile_c_phase_set/publication_test.cpp',
        'tests/support/sw92/test_support.hpp',
        'modules/thermodynamics/include/mpmc/thermodynamics/sw92_phase.hpp',
        'modules/flash/sw92_profile_c_phase_set.md',
        'tests/flash/sw92_profile_c_phase_set/CMakeLists.txt',
        '.github/workflows/sw92_profile_c_phase_set.yml',
    )
    for action in ('opened', 'reopened', 'synchronize', 'ready_for_review'):
        for paths in (['README.md'], ['README.md', 'AGENTS.md', '.github/AGENTS.md']):
            results, ad, thermo = select(paths, action=action)
            assert not any(results.values()) and not ad and not thermo, (
                'root overview edits selected scientific gates', action, paths,
            )
        for path in science_paths:
            expected = select([path], action=action)
            assert expected[0]['sw92_profile_c_phase_set'], ('lost Profile-C owner', path)
            for paths in (['README.md', path], [path, 'README.md']):
                assert select(paths, action=action) == expected, (
                    'root overview changed scientific ownership', action, paths,
                )
        path = 'modules/ad/include/mpmc/ad/math.hpp'
        expected = select([path], action=action)
        assert expected[1], 'AD change must retain its selected suites'
        assert select(['README.md', path], action=action) == expected

def main():
    catalog = json.loads(Path('.github/ci/workflow_map.json').read_text(encoding='utf-8'))
    router_path = '.github/workflows/pr_incremental_ci.yml'
    paths = {p.as_posix() for p in Path('.github/workflows').glob('*.yml')}
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
    # Run the owning MPI CTest directly, clearing cached diagnostic launch flags.
    petsc_run = workflow_run_text('.github/workflows/flow_discretization_petsc.yml')
    assert 'gdb' not in petsc_run
    assert '-DMPIEXEC_PREFLAGS=' in petsc_run
    assert '-U CMAKE_CXX_FLAGS_RELEASE' in petsc_run
    assert petsc_run.count('ctest --test-dir') == 1
    assert '--rerun-failed' not in petsc_run and '|| true' not in petsc_run
    restart_test = Path('tests/flow_discretization/petsc/sw92_transactional_phase_transition_restart_test.cpp').read_text(encoding='utf-8')
    assert 'sw92_transactional_phase_transition_sample6_test(' not in restart_test, 'Sample-6 must execute only from main'
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
        'sw92_transactional_phase_transition_sample6_test.cpp',
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
        r'int\s+main\s*\([^)]*\)\s*\{[\s\S]*?sw92_transactional_phase_transition_sample6_test\s*\(\s*\)\s*;',
        petsc_entry_text,
    ), 'SW92 Sample-6 transactional 2P<->3P regression is compiled but not executed'
    assert re.search(
        r'void\s+headers\s*\(\s*\)\s*\{[\s\S]*?sw92_transactional_phase_transition_restart_header\s*\(\s*\)',
        petsc_entry_text,
    ), 'SW92 transactional restart header probe is compiled but not executed'

    cadence_tests = runpy.run_path('.github/ci/test_cloud_cadence.py')
    cadence_tests['run_tests']()
    # Existing selector regression vectors are run when importing the planner.
    planner = runpy.run_path('.github/ci/plan.py')
    select = planner['select']
    assert_mesh_independent_readers(root, select)
    # Shared 2D geometry affects both importers and their downstream adapters.
    # A deletion is represented by this path too; a rename contributes both sides.
    mesh2d_path = 'modules/mesh/include/mpmc/mesh/linear_cell_mesh_2d.hpp'
    for action in ('opened', 'synchronize', 'ready_for_review'):
        selected, _, _ = select([mesh2d_path], action=action)
        for owner in ('legacy_mesh_core', 'legacy_mesh_external_compatibility', 'legacy_mesh_petsc'):
            assert selected[owner], ('shared 2D builder lost consumer coverage', owner, action)
        assert not selected['legacy_frontend'] and not selected['flow_core']
    before, _, _ = select([mesh2d_path])
    renamed, _, _ = select([mesh2d_path, 'modules/mesh/include/mpmc/mesh/renamed_2d.hpp'])
    assert all(not value or renamed[key] for key, value in before.items())
    unrelated, _, _ = select(['modules/mesh/README.md'])
    assert not unrelated['legacy_mesh_external_compatibility']
    assert_root_readme_selection(select)
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

    assert_root_build_selection(select)

    # Local execution shares the selector and existing test owners. No new CI job.
    local_tests = runpy.run_path('.github/ci/test_local.py')
    local_tests['run_tests']()

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
