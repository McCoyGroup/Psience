
import os

for file in os.listdir(os.path.dirname(os.path.abspath(__file__))):
    if file.endswith("Tests.py"):
        try:
            exec("from .{} import *".format(os.path.splitext(file)[0]) )
        except:
            from traceback import print_exc
            print_exc()