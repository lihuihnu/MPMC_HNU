"""Explicit opt-in dependency preparation for local/CI mesh HDF5 tests.

Pinned upstream source, no vendored code and no system installation. The core
library remains dependency-free when MPMC_MESH_WITH_HDF5 is OFF.
"""
import argparse
import os
from pathlib import Path
import subprocess

REVISION = '7bf340440909d468dbb3cf41f0ea0d87f5050cea'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve(); root.mkdir(parents=True, exist_ok=True)
    source, build, prefix = root/'source', root/'build', root/'install'
    environment = {key.upper(): value for key, value in os.environ.items()} if os.name == 'nt' else os.environ.copy()

    def run(command):
        subprocess.run(command, check=True, env=environment)

    if not source.exists():
        run(['git', 'clone', '--depth', '1', '--branch', 'hdf5_1.14.6',
             'https://github.com/HDFGroup/hdf5.git', str(source)])
    actual = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], env=environment, text=True).strip()
    if actual != REVISION:
        raise RuntimeError(f'HDF5 source revision mismatch: {actual}')
    run(['cmake', '-S', str(source), '-B', str(build), '-DCMAKE_BUILD_TYPE=Release',
         '-DBUILD_TESTING=OFF', '-DBUILD_SHARED_LIBS=OFF', '-DHDF5_BUILD_TOOLS=OFF',
         '-DHDF5_BUILD_EXAMPLES=OFF', '-DHDF5_ENABLE_Z_LIB_SUPPORT=OFF',
         '-DHDF5_ENABLE_SZIP_SUPPORT=OFF', '-DCMAKE_INSTALL_PREFIX='+str(prefix)])
    run(['cmake', '--build', str(build), '--config', 'Release', '--parallel', '2'])
    run(['cmake', '--install', str(build), '--config', 'Release'])
    (root/'revision.txt').write_text(REVISION+'\n', encoding='utf-8')
    print(f'HDF5_ROOT={prefix}')


if __name__ == '__main__':
    main()
