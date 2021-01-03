from Peeves.Doc import *
import os, sys

root = os.path.dirname(os.path.dirname(__file__))
target = os.path.join(root, "docs")
sys.path.insert(0, root)
doc_config = {
    "config": {
        "title": "Psience Dev Branch Documentation",
        "path": "Psience",
        "url": "https://mccoygroup.github.io/Psience/",
        "gh_username": "McCoyGroup",
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
    "readme": os.path.join(root, "README.md")
}
DocBuilder(**doc_config).build()