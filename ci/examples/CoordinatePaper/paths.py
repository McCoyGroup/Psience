

import os

root_dir = os.path.expanduser("~/Documents/Postdoc/Projects/CoordinatePaper")
def path(*bits):
    return os.path.join(root_dir, *bits)
def torsion_scan_path(*bits):
    return path('torsion_scan', *bits)