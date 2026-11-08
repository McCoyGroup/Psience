"""Test discovery is explicit; importing this package does not import tests."""


def load_tests(loader, tests, pattern):
    import os
    import sys

    # Peeves selects test modules before importing them, just as it selects
    # methods afterwards through each module's own load_tests hook.
    harness = sys.modules.get("Peeves.TestUtils")
    manager = getattr(getattr(harness, "ManagedTestLoader", None), "manager", None)
    selected = None if manager is None else manager.test_files
    if selected == "All":
        selected = None
    elif isinstance(selected, str):
        selected = selected.split(",")

    for filename in sorted(os.listdir(os.path.dirname(__file__))):
        if not filename.endswith("Tests.py"):
            continue
        name = filename[:-3]
        if selected is not None and name not in selected:
            continue
        qualified_name = __name__ + "." + name
        module_tests = loader.loadTestsFromName(qualified_name)
        module = sys.modules.get(qualified_name)
        if manager is not None and module is not None and not hasattr(module, "load_tests"):
            # Plain unittest modules need the same Peeves queueing as modules
            # which import Peeves' load_tests hook themselves.
            module_tests = harness.load_tests(loader, module_tests, pattern)
        tests.addTests(module_tests)
    return tests
