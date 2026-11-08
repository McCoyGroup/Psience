"""Test collection and core checks must not load optional integrations."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import textwrap
import unittest

PACKAGE = 'Psience'
ROOT = Path(__file__).resolve().parents[2]


class TestDependencyIsolationTests(unittest.TestCase):
    def run_child(self, code, extras=False):
        environment = dict(os.environ, PYTHONDONTWRITEBYTECODE="1",
                           PYTHONNOUSERSITE="1",
                           PYTHONPATH=os.pathsep.join([str(ROOT)] + [os.path.abspath(p) for p in sys.path]))
        environment.pop("MCUTILS_TEST_EXTRAS", None)
        if extras:
            environment["MCUTILS_TEST_EXTRAS"] = "1"
        with tempfile.TemporaryDirectory() as directory:
            process = subprocess.run([sys.executable, "-B", "-c", textwrap.dedent(code)],
                                     cwd=directory, env=environment, capture_output=True,
                                     text=True, timeout=60)
        self.assertEqual(process.returncode, 0, process.stdout + process.stderr)

    def test_importing_test_package_does_not_collect_tests(self):
        self.run_child(f'''
            import importlib.util
            import sys
            spec = importlib.util.spec_from_file_location(
                "_isolation_tests", {str(ROOT / "ci/tests/__init__.py")!r},
                submodule_search_locations=[{str(ROOT / "ci/tests")!r}])
            module = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = module
            spec.loader.exec_module(module)
            assert not any(name.startswith(spec.name + ".") for name in sys.modules)
            assert callable(module.load_tests)
        ''')

    def test_selected_discovery_does_not_import_unrelated_modules(self):
        self.run_child(f'''
            import importlib.util
            import pathlib
            import shutil
            import sys
            import types
            import unittest
            directory = pathlib.Path.cwd()
            shutil.copyfile({str(ROOT / "ci/tests/__init__.py")!r}, directory / "__init__.py")
            (directory / "CoreTests.py").write_text(
                "import unittest\\nclass CoreTests(unittest.TestCase):\\n"
                "    def test_core(self): self.assertEqual(2 + 2, 4)\\n")
            (directory / "UnrelatedTests.py").write_text(
                "import builtins\\nbuiltins.unrelated_test_imported = True\\n"
                "raise AssertionError('unrelated test imported')\\n")
            harness = types.ModuleType("Peeves.TestUtils")
            harness.queued = []
            harness.load_tests = lambda loader, tests, pattern: (harness.queued.append(tests.countTestCases()) or tests)
            harness.ManagedTestLoader = types.SimpleNamespace(
                manager=types.SimpleNamespace(test_files="CoreTests"))
            sys.modules["Peeves.TestUtils"] = harness
            spec = importlib.util.spec_from_file_location(
                "_isolation_tests", directory / "__init__.py",
                submodule_search_locations=[str(directory)])
            module = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = module
            spec.loader.exec_module(module)
            suite = unittest.TestLoader().loadTestsFromModule(module)
            assert suite.countTestCases() == 1
            assert harness.queued == [1], harness.queued
            result = unittest.TestResult()
            suite.run(result)
            assert result.wasSuccessful(), (result.errors, result.failures)
            import builtins
            assert not hasattr(builtins, "unrelated_test_imported")
        ''')

    def test_optional_tests_are_skipped_and_core_checks_still_run(self):
        target = "ColorPaletteTests.py" if PACKAGE == "McUtils" else "PolarizabilityModeTests.py"
        self.run_child(f'''
            import importlib.abc
            import importlib.util
            import sys
            import unittest
            optional = {{"torch", "markdown", "sympy", "matplotlib", "ipywidgets",
                        "IPython", "h5py", "rdkit", "pyscf", "ase", "aimnet",
                        "mace", "plotly", "pyinstrument", "nglview"}}
            requested = []
            class NoOptional(importlib.abc.MetaPathFinder):
                def find_spec(self, fullname, path=None, target=None):
                    if fullname.split(".")[0] in optional:
                        requested.append(fullname)
                        raise ModuleNotFoundError("Optional dependency requested: " + fullname,
                                                  name=fullname)
            sys.meta_path.insert(0, NoOptional())
            spec = importlib.util.spec_from_file_location(
                "_isolated_module", {str(ROOT / "ci/tests" / target)!r})
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            result = unittest.TestResult()
            unittest.TestLoader().loadTestsFromModule(module).run(result)
            assert result.wasSuccessful(), (result.errors, result.failures)
            assert result.testsRun > len(result.skipped), result.testsRun
            assert result.skipped
            assert not requested, requested
            assert not optional.intersection(name.split(".")[0] for name in sys.modules)
        ''')

    def test_explicit_opt_in_does_not_eagerly_import_backend(self):
        target = "ClaudePlotsTests.py" if PACKAGE == "McUtils" else "PolarizabilityModeTests.py"
        class_name = "ClaudePlotsTests" if PACKAGE == "McUtils" else "AIMNetResponseModeTests"
        backend = "matplotlib" if PACKAGE == "McUtils" else "torch"
        self.run_child(f'''
            import importlib.util
            import sys
            spec = importlib.util.spec_from_file_location(
                "_isolated_module", {str(ROOT / "ci/tests" / target)!r})
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            assert not getattr(module.{class_name}, "__unittest_skip__", False)
            assert {backend!r} not in sys.modules
        ''', extras=True)


if __name__ == "__main__":
    unittest.main()
