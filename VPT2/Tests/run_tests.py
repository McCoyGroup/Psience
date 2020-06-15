"""Misc tests for my misc packages
Each package can package its own tests I figure?
"""

import os, sys

test_dir = os.path.dirname(os.path.abspath(__file__))
test_root = os.path.dirname(test_dir)
sys.path.insert(0, test_root)

from Peeves.TestUtils import TestManager
TestManager.base_dir = test_root

# provide a nice way to automatically pipe print output to stderr so it appears in the regular
# output area for the unit tests
stdout = sys.stdout
try:
    sys.stdout = sys.stderr
    TestManager.run()
finally:
    sys.stdout = stdout