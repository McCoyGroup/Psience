"""
Psience
A simple package for working with quantum chemistry problems
"""

# pulled from Ryan's stuff, lightly modified
from setuptools import setup, find_packages

short_description = __doc__.split("\n")

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except:
    long_description = "\n".join(short_description[2:])

def get_version():
    import subprocess
    run_out = subprocess.run(['git', 'describe', '--tags', '--abbrev=0'], capture_output=True)
    return run_out.stdout.decode().strip().strip("v")

setup(
    # Self-descriptive entries which should always be present
    name='mccoygroup-psience',
    author='Mark Boyer',
    author_email='b3m2a1@uw.edu',
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=get_version(),
    license='MIT',

    # Which Python importable modules should be included when your package is installed
    # Handled automatically by setuptools. Use 'exclude' to prevent some specific
    # subpackage(s) from being added, if needed
    packages=find_packages(),

    # Optional include package data to ship with your package
    # Customize MANIFEST.in if the general case does not suit your needs
    # Comment out this line to prevent the files from being packaged with your software
    # include_package_data=True
    install_requires=[
        'mccoygroup-mcutils>=0.0.9.3'
    ]
)