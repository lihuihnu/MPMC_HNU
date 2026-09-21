from pathlib import Path
p = Path('/tmp/entry_migration.py')
s = p.read_text()
def replace(a, b):
    global s
    assert a in s, a
    s = s.replace(a, b)
replace("'pull_request', 'workflow_call', 'workflow_dispatch'}", "'pull_request', 'push', 'workflow_call', 'workflow_dispatch'}")
replace("if 'pull_request' not in events:", "if not ({'pull_request', 'push'} & set(events)):\n")
replace("wf['on'].pop('pull_request')", "wf['on'].pop('pull_request', None)\n    wf['on'].pop('push', None)")
replace("s = s.replace('github.workflow',", "s = re.sub(r'github\\.workflow(?![A-Za-z0-9_])',")
replace("original.get('name', path).replace(\"'\", \"''\") + \"'\")", "original.get('name', path).replace(\"'\", \"''\") + \"'\", s)")
replace("github.event_name == 'pull_request' && needs.impact.outputs.trusted", "(github.event_name == 'pull_request' || github.event_name == 'push') && needs.impact.outputs.trusted")
replace("selector = spec['selector']\n        eligible", "selector = spec['selector']\n        if 'pull_request' not in spec['original_events']:\n            result[key] = False\n            continue\n        eligible")
replace("    else:\n        result = {k: True for k in RULES['ownership']}\n", "    elif os.environ.get('GITHUB_EVENT_NAME') == 'push':\n        head = event['after']\n        before = event.get('before', '')\n        if re.fullmatch(r'[0-9a-f]{40}', before) and before != '0' * 40 and ancestor(before, head):\n            raw = subprocess.check_output(['git', 'diff', '--name-only', '--no-renames', '-z', before, head])\n        else:\n            raw = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', '-z', head])\n        paths = [p for p in raw.decode('utf-8').split('\\0') if p]\n        result = {k: False for k in RULES['ownership']}\n        ad, thermo = [], []\n        branch = event.get('ref', '').removeprefix('refs/heads/')\n        for spec in MAP['workflows'].values():\n            if 'route_key' not in spec:\n                continue\n            selector = spec['original_events'].get('push')\n            eligible = selector is not None\n            selector = selector or {}\n            eligible = eligible and matches(branch, selector.get('branches', ['**'])) and not matches(branch, selector.get('branches-ignore', []))\n            result[spec['route_key']] = eligible and any(matches(p, selector.get('paths', ['**'])) and not matches(p, selector.get('paths-ignore', [])) for p in paths)\n    else:\n        result = {k: True for k in RULES['ownership']}\n")
replace("root['permissions'] = {'contents': 'read', 'actions': 'read'}", "push_specs = [sp['original_events']['push'] or {} for sp in registry['workflows'].values() if 'push' in sp['original_events']]\nif push_specs:\n    assert all(set(sp) <= {'branches', 'branches-ignore', 'paths', 'paths-ignore'} for sp in push_specs), push_specs\n    root['on']['push'] = {'branches': sorted({b for sp in push_specs for b in sp.get('branches', ['**'])}), 'paths': sorted({q for sp in push_specs for q in sp.get('paths', ['**'])})}\nprint('PRESERVED_PUSH_SELECTORS', json.dumps(push_specs))\nroot['permissions'] = {'contents': 'read', 'actions': 'read'}")
# Dump all original entry kinds before the exact semantic guards run.
replace("root = copy.deepcopy(workflows[ROUTER])", "print('ORIGINAL_EVENT_INVENTORY', json.dumps({p: w.get('on') for p,w in workflows.items()}, ensure_ascii=False))\nroot = copy.deepcopy(workflows[ROUTER])")
p.write_text(s)
