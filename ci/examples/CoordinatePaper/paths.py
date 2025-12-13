

import os

root_dir = os.path.expanduser("~/Documents/Postdoc/Projects/CoordinatePaper")
def path(*bits):
    return os.path.join(root_dir, *bits)
torsion_scan_root = 'torsion_scan'
def torsion_scan_path(*bits):
    return path(torsion_scan_root, *bits)