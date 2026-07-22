
import os, sys
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
target = os.path.join(root, "ci", "stubs")
sys.path.insert(0, root)

if __name__ == "__main__":
    from McUtils.Docs import *
    import Psience
    StubSummaryBuilder(verbose=True, out_dir=target, allow_static_mode=False).generate_all(Psience)