# Windows Electron desktop installer gate

This gate packages the already-audited Windows Electron preview and its embedded
native PT staging tree into one WiX MSI. It does not add EOS, flash, fallback,
retry, or result-interpretation logic.

The CI acceptance path is deliberately stronger than administrative MSI
extraction:

1. rebuild the exact Windows x64 native staging payload from the locked Conan
   graph with `--no-remote --build=never`;
2. build the shared React/Electron renderer and package the existing desktop
   preview with `@electron/packager`;
3. wrap that complete directory in a per-machine WiX MSI under
   `Program Files\MPMC PT Desktop`;
4. perform a real silent `msiexec /i` installation;
5. launch the **installed** Electron executable in an opt-in hidden smoke mode;
6. feed a versioned external smoke plan through the existing preload bridge and
   require one scientific result from each configured PR76, SW92 and CPA
   backend;
7. verify the Start Menu shortcut resolves to the installed executable; and
8. perform `msiexec /x` and require both the executable and shortcut to be
   removed.

The installed-smoke plan is supplied only by CI. The production Electron bundle
contains generic plan parsing and bridge orchestration but no model-specific
backend IDs, compositions, EOS branching, or scientific acceptance thresholds.

The MSI and embedded binaries remain unsigned engineering candidates. A green
gate establishes installability and functional desktop integration on the
official Windows runner; it does not establish Authenticode trust or public
release eligibility.
