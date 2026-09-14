from __future__ import annotations

import pathlib
import unittest

INSTALLER_ROOT = pathlib.Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = INSTALLER_ROOT.parents[2]


class CrossPlatformPreviewPackagingContractTest(unittest.TestCase):
    def test_preview_executable_suffix_is_platform_specific(self):
        packager = (
            REPOSITORY_ROOT / "frontend" / "desktop" / "packageDesktop.mjs"
        ).read_text(encoding="utf-8")
        self.assertIn(
            "platform === 'win32' ? `${electronExecutable}.exe` : electronExecutable",
            packager,
        )
        self.assertNotIn(
            "releaseIdentity?.releaseExecutable ?? 'MPMC-PT-Desktop-Preview.exe'",
            packager,
        )
        self.assertIn("releaseIdentity?.releaseExecutable", packager)


if __name__ == "__main__":
    unittest.main()
