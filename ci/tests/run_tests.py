"""
Tests for things in the Psience packages
"""

import os, sys
run_file = os.path.abspath(__file__)
test_dir = os.path.dirname(run_file)
test_root = os.path.dirname(test_dir)
test_pkg = os.path.basename(test_dir)
sys.path.insert(0, os.path.dirname(test_root))

from Peeves.TestUtils import TestManager
TestManager.test_root = test_root
TestManager.test_pkg = test_pkg

# provide a nice way to automatically pipe print output to stderr so it appears in the regular
# output area for the unit tests
stdout = sys.stdout
try:
    sys.stdout = sys.stderr
    TestManager.run()
finally:
    sys.stdout = stdout