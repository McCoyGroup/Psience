"""These are the set of tests to the PyVPT package
Every type of tests should be its own module and should be tagged as either a debugTest, validationTest, or timingTest


"""

if __name__ == "__main__":
    # we'll put the parent dir on path
    import sys, os
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    sys.path.insert(0, base_dir)
    from PyVPT.Peeves import TestManager
    TestManager.run()
