from __future__ import annotations

import pathlib
import unittest

INSTALLER_ROOT = pathlib.Path(__file__).resolve().parents[1]


class WindowsGuidCMakeCompatibilityTest(unittest.TestCase):
    def test_guid_validation_avoids_unsupported_counted_repetition(self):
        cmake = (INSTALLER_ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
        self.assertNotIn("[0-9A-F]{8}", cmake)
        self.assertNotIn("[0-9A-F]{4}", cmake)
        self.assertNotIn("[0-9A-F]{12}", cmake)
        self.assertIn("CMake's regular-expression engine does not support {n} repetition", cmake)
        self.assertIn('NOT mpmc_upgrade_guid MATCHES "${mpmc_guid_regex}"', cmake)
        self.assertIn('NOT mpmc_rc_product_guid MATCHES "${mpmc_guid_regex}"', cmake)


if __name__ == "__main__":
    unittest.main()
