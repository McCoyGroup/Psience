from McUtils.Docs import *
import os, sys

root = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
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
            "id": "Psience",
            'tests_root': os.path.join(root, "ci", "tests")
        }
    ],
    "root": root,
    "target": target,
    "readme": os.path.join(root, "blurb.md"),
    'templates_directory': os.path.join(root, 'ci', 'docs', 'templates'),
    'examples_directory': os.path.join(root, 'ci',  'docs', 'examples')
}
ExamplesParser.IGNORE_UNHANDLED_STATEMENTS = True
DocBuilder(**doc_config).build()