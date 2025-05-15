from McUtils.Docs import *
import os, sys

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
target = os.path.join(root, "ci", "docs")
sys.path.insert(0, root)
doc_config = {
    "config": {
        "title": "Psience Documentation",
        "path": "Psience",
        "url": "https://mccoygroup.github.io/Psience/",
        "gh_username": "McCoyGroup",
        "gh_repo": "Psience",
        "gh_branch": "master",
        "footer": "Brought to you by the McCoy Group"
    },
    "packages": [
        {
            "id": "Psience"
        }
    ],
    "root": root,
    "target": target,
    "readme": os.path.join(root, "blurb.md"),
    'templates_directory': os.path.join(root, 'ci', 'docs', 'templates'),
    'examples_directory': os.path.join(root, 'ci',  'docs', 'examples'),
    'tests_directory': os.path.join(root, "ci", "tests")
}

ExamplesParser.IGNORE_UNHANDLED_STATEMENTS = True
DocBuilder(**doc_config).build()