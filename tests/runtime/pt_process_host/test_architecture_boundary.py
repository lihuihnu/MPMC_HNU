import pathlib
import unittest


ROOT = pathlib.Path(__file__).parents[3]


class ArchitectureBoundaryTest(unittest.TestCase):
    def test_host_and_adapter_have_no_model_specific_dependencies(self):
        files = (
            ROOT / "modules/pt_process/src/process_host.cpp",
            ROOT / "modules/pt_process/include/mpmc/pt_process/process_host.hpp",
            ROOT / "modules/pt_process/src/composition_root.cpp",
            ROOT / "modules/pt_process/include/mpmc/pt_process/composition_root.hpp",
            ROOT / "modules/runtime_grpc/src/pt_grpc_adapter.cpp",
            ROOT / "modules/runtime_grpc/include/mpmc/runtime_grpc/pt_grpc_adapter.hpp",
        )
        for path in files:
            text = path.read_text(encoding="utf-8").lower()
            with self.subTest(path=path):
                for model_name in ("pr76", "sw92", "cpa"):
                    self.assertNotIn(model_name, text)

    def test_adapter_has_one_scientific_delegation(self):
        adapter = (
            ROOT / "modules/runtime_grpc/src/pt_grpc_adapter.cpp"
        ).read_text(encoding="utf-8")
        self.assertEqual(adapter.count("service_.solve(service_request)"), 1)

    def test_typed_factories_construct_but_never_solve(self):
        factory = (
            ROOT / "modules/pt_process/src/configured_backends.cpp"
        ).read_text(encoding="utf-8")
        self.assertNotIn(".solve(", factory)
        self.assertNotIn("solve_pr76", factory)
        self.assertNotIn("solve_sw92", factory)
        self.assertNotIn("solve_cpa", factory)


if __name__ == "__main__":
    unittest.main()
