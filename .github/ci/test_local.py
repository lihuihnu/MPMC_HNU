"""Local planner/executor regressions, owned by the existing CI impact job."""
from contextlib import redirect_stdout
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('mpmc_local', Path(__file__).with_name('local.py'))
local = importlib.util.module_from_spec(spec)
spec.loader.exec_module(local)


class LocalTests(unittest.TestCase):
    def test_documentation_does_not_select_science(self):
        plan = local.make_plan(['README.md', 'AGENTS.md', 'docs/development.md'])
        self.assertEqual(plan['local_suites'], ['governance'])
        self.assertEqual(plan['cloud_only_gates'], [])

    def test_local_tools_owned_by_governance(self):
        for path in ['.github/ci/local.py', '.github/ci/test_local.py', '.github/ci/requirements-local.txt']:
            self.assertEqual(local.make_plan([path])['local_suites'], ['governance'])

    def test_windows_paths_and_unknown_owner(self):
        self.assertEqual(local.make_plan([r'tests\ad\math\math_test.cpp']),
                         local.make_plan(['tests/ad/math/math_test.cpp']))
        with self.assertRaisesRegex(ValueError, 'lacks cloud ownership'):
            local.make_plan(['tests/new_owner/new_test.cpp'])
        for path in ['../outside.cpp', 'C:/outside.cpp', '/outside.cpp']:
            with self.assertRaises(ValueError):
                local.make_plan([path])

    def test_downstream_and_unsupported_gates_are_preserved(self):
        plan = local.make_plan(['modules/mesh/include/mpmc/mesh/topology.hpp'])
        self.assertIn('legacy_mesh_core', plan['local_suites'])
        self.assertIn('legacy_mesh_petsc', plan['cloud_only_gates'])
        self.assertIn('ad.arithmetic', plan['local_suites'])

    def test_router_semantic_edits_require_cloud_validation(self):
        plan = local.make_plan(['.github/workflows/pr_incremental_ci.yml'])
        self.assertIn('central_router_semantic_validation', plan['cloud_only_gates'])

    def test_recipes_use_current_workflows(self):
        for suite in local.NATIVE:
            source, targets, pattern, sanitizer = local.native_recipe(suite)
            self.assertTrue((local.ROOT / source / 'CMakeLists.txt').exists())
            self.assertTrue(targets and pattern and sanitizer)
        with patch.object(local, 'workflow_text', return_value='changed layout'):
            with self.assertRaisesRegex(ValueError, 'requires review'):
                local.native_recipe('flow_core')

    def test_ad_arithmetic_and_thermo_keep_extra_contracts(self):
        root = Path(tempfile.gettempdir()) / 'mpmc-command-preview'
        commands = local.commands('ad.arithmetic', root, 'Release', 2)
        self.assertTrue(any('root_libraries/verify.py' in ' '.join(c).replace('\\', '/') for c in commands))
        self.assertTrue(any('tests/ad/standalone' in ' '.join(c).replace('\\', '/') for c in commands))
        commands = local.commands('thermo.pr76_pt', root, 'Release', 2)
        self.assertTrue(any('reference_decimal.py' in ' '.join(c) for c in commands))
        for suite in local.SUITES:
            commands = local.commands(suite, root, 'Release', 2)
            for command in commands:
                if command[0] == 'cmake' and '-B' in command:
                    folder = Path(command[command.index('-B') + 1])
                    self.assertLessEqual(len(folder.name), 23)
                if command[0] == 'ctest':
                    self.assertIn('--no-tests=error', command)

    def test_output_must_be_outside_checkout(self):
        with self.assertRaises(ValueError):
            local.output_root(local.ROOT / 'build')

    def test_cumulative_diff_staged_unstaged_rename_delete_untracked(self):
        with tempfile.TemporaryDirectory(prefix='mpmc-local-git-') as directory:
            root = Path(directory)
            def git(*args):
                return subprocess.check_output(['git', '-c', 'user.name=Local test',
                    '-c', 'user.email=local-test@example.invalid', *args], cwd=root, stderr=subprocess.DEVNULL)
            git('init', '-b', 'main')
            for name in ('old name.cpp', 'deleted.cpp', 'unstaged.cpp'):
                (root / name).write_text('original\n', encoding='utf-8')
            git('add', '.')
            git('commit', '-m', 'baseline')
            git('branch', 'baseline')
            git('mv', 'old name.cpp', 'new name.cpp')
            git('rm', 'deleted.cpp')
            git('commit', '-m', 'committed change')
            (root / 'staged.cpp').write_text('staged\n', encoding='utf-8')
            git('add', 'staged.cpp')
            (root / 'unstaged.cpp').write_text('modified\n', encoding='utf-8')
            (root / '未跟踪.cpp').write_text('untracked\n', encoding='utf-8')
            with patch.object(local, 'ROOT', root):
                _, paths = local.changed_paths('baseline')
            self.assertEqual(set(paths), {'old name.cpp', 'new name.cpp', 'deleted.cpp',
                                         'staged.cpp', 'unstaged.cpp', '未跟踪.cpp'})

    def execute_stub(self, exit_code, cloud_only=(), changed=False):
        with tempfile.TemporaryDirectory(prefix='mpmc-local-result-') as directory:
            plan = {'local_suites': ['governance'], 'cloud_only_gates': list(cloud_only)}
            snapshots = [{'head': 'before'}, {'head': 'after' if changed else 'before'}]
            with patch.object(local, 'snapshot', side_effect=snapshots), \
                 patch.object(local, 'doctor', return_value={}), \
                 patch.object(local, 'commands', return_value=[[sys.executable, '-c', f'raise SystemExit({exit_code})']]), \
                 redirect_stdout(io.StringIO()):
                code = local.execute(plan, Path(directory), 'Release', 2)
            reports = list(Path(directory).glob('runs/*/result.json'))
            self.assertEqual(len(reports), 1)
            return code, json.loads(reports[0].read_text(encoding='utf-8'))

    def test_failure_is_not_success(self):
        code, report = self.execute_stub(7)
        self.assertEqual(code, 1)
        self.assertEqual(report['status'], 'failed')
        self.assertEqual(report['results'][0]['exit_code'], 7)

    def test_cloud_only_means_partial(self):
        code, report = self.execute_stub(0, ['flow_discretization_petsc'])
        self.assertEqual(code, 2)
        self.assertEqual(report['status'], 'partial')

    def test_source_change_invalidates_result(self):
        code, report = self.execute_stub(0, changed=True)
        self.assertEqual(code, 1)
        self.assertIn('Source changed', report['error'])


def run_tests():
    result = unittest.TextTestRunner(verbosity=1).run(unittest.defaultTestLoader.loadTestsFromTestCase(LocalTests))
    if not result.wasSuccessful():
        raise AssertionError('Local developer workflow regressions failed')


if __name__ == '__main__':
    run_tests()
