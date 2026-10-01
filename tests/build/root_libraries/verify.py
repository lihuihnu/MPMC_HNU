"""Root build contract; owned once per platform by the AD arithmetic CI job."""
import argparse
from itertools import product
import json
from pathlib import Path
import subprocess
import tempfile


def run(*args):
    subprocess.run(args, check=True)


def tests_in(build, config):
    result = subprocess.check_output(
        ['ctest', '--test-dir', str(build), '-C', config, '--show-only=json-v1'],
        text=True,
    )
    return sorted(test['name'] for test in json.loads(result)['tests'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default='Debug')
    parser.add_argument('--compiler', default='')
    parser.add_argument('--sanitizer', choices=['ON', 'OFF'], default='OFF')
    args = parser.parse_args()
    source = Path(__file__).resolve().parent
    root = source.parents[2]
    common = [f'-DCMAKE_BUILD_TYPE={args.config}']
    if args.compiler:
        common.append(f'-DCMAKE_CXX_COMPILER={args.compiler}')
    options = ('MPMC_ENABLE_THERMODYNAMICS', 'MPMC_ENABLE_FLASH',
               'MPMC_ENABLE_MESH', 'MPMC_ENABLE_DISCRETIZATION',
               'MPMC_ENABLE_FLOW_DISCRETIZATION', 'MPMC_ENABLE_WELL_DISCRETIZATION')
    with tempfile.TemporaryDirectory(prefix='mpmc-root-build-') as temporary:
        work = Path(temporary)
        # Root defaults must still register exactly the original AD test.
        default_build = work / 'default'
        run('cmake', '-S', str(root), '-B', str(default_build), *common)
        assert tests_in(default_build, args.config) == ['ad.dual']
        cache = (default_build / 'CMakeCache.txt').read_text(encoding='utf-8')
        for option in options:
            assert f'{option}:BOOL=OFF' in cache, (option, 'default changed')
        # Exercise the documented preset while keeping outputs outside the repo.
        preset_build = work / 'preset'
        subprocess.run(
            ['cmake', '--preset', 'flash-libraries', '-B', str(preset_build), *common],
            cwd=root, check=True,
        )
        preset_cache = (preset_build / 'CMakeCache.txt').read_text(encoding='utf-8')
        for setting in ('BUILD_TESTING:BOOL=OFF', 'MPMC_ENABLE_FLASH:BOOL=ON',
                        'MPMC_ENABLE_MESH:BOOL=OFF', 'MPMC_ENABLE_DISCRETIZATION:BOOL=OFF',
                        'MPMC_ENABLE_FLOW_DISCRETIZATION:BOOL=OFF',
                        'MPMC_ENABLE_WELL_DISCRETIZATION:BOOL=OFF'):
            assert setting in preset_cache, ('library preset drifted', setting)
        assert tests_in(preset_build, args.config) == []
        run('cmake', '--build', str(preset_build), '--config', args.config)
        # Exhaust the original chains with/without the well bridge, including
        # explicitly disabled lower modules required by the enabled bridge.
        for values in product(('OFF', 'ON'), repeat=len(options)):
            build = work / '-'.join(values)
            run('cmake', '-S', str(source), '-B', str(build), *common,
                *(f'-D{option}={value}' for option, value in zip(options, values)),
                f'-DMPMC_AD_ENABLE_SANITIZERS={args.sanitizer}')
            cache = (build / 'CMakeCache.txt').read_text(encoding='utf-8')
            for option, value in zip(options, values):
                assert f'{option}:BOOL={value}' in cache, (option, 'caller option changed')
            expected = ['build.root_libraries']
            if 'ON' in values[2:4]:
                expected.append('build.root_spatial_libraries')
            if 'ON' in values[4:6]:
                expected.extend(['build.root_flow_thermodynamics', 'build.root_flow_discretization'])
            if values[5] == 'ON':
                expected.extend(['build.root_well', 'build.root_well_discretization'])
            assert tests_in(build, args.config) == sorted(expected)
            run('cmake', '--build', str(build), '--config', args.config, '--parallel', '2')
            run('ctest', '--test-dir', str(build), '-C', args.config,
                '--output-on-failure', '--no-tests=error')
        # Disabling modules in an existing build must remove their targets, too.
        build = work / '-'.join('ON' for _ in options)
        run('cmake', '-S', str(source), '-B', str(build), *common,
            *(f'-D{option}=OFF' for option in options))
        assert tests_in(build, args.config) == ['build.root_libraries']
        run('cmake', '--build', str(build), '--config', args.config, '--parallel', '2')
        run('ctest', '--test-dir', str(build), '-C', args.config,
            '--output-on-failure', '--no-tests=error')
    print('ROOT_BUILD_CONTRACT_OK defaults; library preset; 64 option combinations; reconfigure')


if __name__ == '__main__':
    main()
