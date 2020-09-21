from distutils.core import setup, Extension
import shutil, os, sys

lib_dir = os.path.dirname(os.path.abspath(__file__))
libname="constantPot"

def get_extension():
    return Extension(
        libname,
        sources = [ libname+'.cpp' ]
    )

def setup_compile():
    return setup(
        name = libname,
        version = '1.0',
        description = 'A C potential to hook into python',
        ext_modules = [get_extension()],
        language = "compiled++"
    )

def compile():

    curdir = os.getcwd()
    try:
        os.chdir(lib_dir)
        sysargv1 = sys.argv
        sys.argv = ['build', 'build_ext', '--inplace']
        try:
            setup_compile()
        finally:
            sys.argv = sysargv1
    finally:
        os.chdir(curdir)

def find_source():
    ext = ""
    target = os.path.join(lib_dir, libname)
    src = None
    for f in os.listdir(lib_dir):
        if f.startswith(libname) and f.endswith(".so"):
            ext = ".so"
            src = os.path.join(lib_dir, f)
            target += ext
            break
        elif f.startswith(libname) and f.endswith(".pyd"):
            ext = ".pyd"
            src = os.path.join(lib_dir, f)
            target += ext
            break

    return src, target

def load():
    compile()
    src, target = find_source()
    if src is not None:
        # try:
        #     os.remove(target)
        # except:
        #     pass
        # os.rename(src, target)
        shutil.rmtree(os.path.join(lib_dir, "build"))
        return src
    else:
        return None

if __name__=="__main__":
    loaded = load()
    if loaded is None:
        raise ImportError("{}: library didn't compile".format(libname))