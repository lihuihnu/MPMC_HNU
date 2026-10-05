"""Local-first cadence regressions, owned by verify_workflows.py (no new job)."""
from contextlib import redirect_stdout
import importlib.util
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import yaml

spec = importlib.util.spec_from_file_location('mpmc_cadence_plan', Path(__file__).with_name('plan.py'))
planner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(planner)


class CloudCadenceTests(unittest.TestCase):
    def test_draft_gate_and_dependency_closure(self):
        workflow = yaml.safe_load(Path('.github/workflows/pr_incremental_ci.yml').read_text(encoding='utf-8'))
        jobs = workflow['jobs']
        gate = "github.event_name != 'pull_request' || github.event.pull_request.draft == false"
        self.assertEqual(jobs['impact']['if'], '${{ ' + gate + ' }}')
        self.assertEqual(jobs['result']['if'], '${{ always() && (' + gate + ') }}')
        for name, job in jobs.items():
            if name in ('impact', 'result'):
                continue
            with self.subTest(job=name):
                needs = job.get('needs', [])
                self.assertIn('impact', [needs] if isinstance(needs, str) else needs)
                self.assertIn("needs.impact.outputs.trusted == 'true'", job['if'])
                # Legacy jobs may use always(), but it must remain behind the
                # false trusted output when impact was skipped on a draft.
                self.assertRegex(job['if'], r"^(?:\(github.event_name == 'pull_request' \|\| github.event_name == 'push'\) && )?needs\.impact\.outputs\.trusted == 'true' && ")
        events = workflow.get('on', workflow.get(True))
        self.assertIn('ready_for_review', events['pull_request']['types'])
        self.assertIn('converted_to_draft', events['pull_request']['types'])
        # Returning to draft supersedes an in-flight formal run on the same PR.
        self.assertTrue(workflow['concurrency']['cancel-in-progress'])
        self.assertIn('github.ref', workflow['concurrency']['group'])
        self.assertIn('workflow_dispatch', events)
        self.assertIn('push', events)

    def run_pr(self, action, draft, trusted=True):
        repository = 'example/project'
        event = dict(action=action, pull_request=dict(
            draft=draft, number=7, base=dict(sha='base', ref='main'),
            head=dict(sha='head', ref='feature', repo=dict(full_name=repository if trusted else 'fork/project'))))
        commands = []

        def git_output(command, **_):
            commands.append(command)
            if command == ['git', 'merge-base', 'base', 'head']:
                return 'merge-base\n'
            if command == ['git', 'diff', '--name-only', '--no-renames', '-z', 'merge-base', 'head']:
                # Represents an early feature slice, even if recent commits only change docs.
                return b'modules/mesh/include/mpmc/mesh/linear_cell_mesh_2d.hpp\0'
            self.fail(f'unexpected command: {command}')

        with tempfile.TemporaryDirectory(prefix='mpmc-cadence-') as directory:
            directory = Path(directory)
            source, output = directory/'event.json', directory/'output.txt'
            source.write_text(json.dumps(event), encoding='utf-8')
            env = dict(GITHUB_EVENT_PATH=str(source), GITHUB_OUTPUT=str(output),
                       GITHUB_EVENT_NAME='pull_request', GITHUB_REPOSITORY=repository)
            stream = io.StringIO()
            with patch.dict(os.environ, env), redirect_stdout(stream), \
                    patch.object(planner.subprocess, 'check_output', side_effect=git_output), \
                    patch.object(planner.urllib.request, 'urlopen', side_effect=AssertionError('PR must not query checkpoints')):
                if draft is not False:
                    with self.assertRaisesRegex(RuntimeError, 'Draft PR'):
                        planner.main()
                    self.assertFalse(output.exists())
                    self.assertEqual(commands, [])
                    return
                planner.main()
            self.assertIn('cumulative PR diff (feature acceptance)', stream.getvalue())
            values = dict(line.split('=', 1) for line in output.read_text(encoding='utf-8').splitlines())
            self.assertEqual(values['trusted'], str(trusted).lower())
            for gate in ('legacy_mesh_core', 'legacy_mesh_external_compatibility', 'legacy_mesh_petsc'):
                self.assertEqual(values[gate], str(trusted).lower())

    def test_draft_and_missing_state_fail_closed(self):
        for action in ('opened', 'synchronize', 'converted_to_draft'):
            for draft in (True, None):
                with self.subTest(action=action, draft=draft):
                    self.run_pr(action, draft)

    def test_ready_and_acceptance_repairs_cover_entire_feature(self):
        for action in ('opened', 'reopened', 'ready_for_review', 'synchronize'):
            with self.subTest(action=action):
                self.run_pr(action, False)

    def test_untrusted_fork_still_cannot_schedule_private_jobs(self):
        self.run_pr('ready_for_review', False, trusted=False)


def run_tests():
    result = unittest.TextTestRunner(verbosity=1).run(unittest.defaultTestLoader.loadTestsFromTestCase(CloudCadenceTests))
    if not result.wasSuccessful():
        raise AssertionError('cloud cadence regressions failed')


if __name__ == '__main__':
    run_tests()
