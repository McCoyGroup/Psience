#!/bin/bash
set -e

##
## Runs CI
##

# get into the parent folder and merge in changes from the master branch
cd /home/Psience
git config user.name ${GITHUB_ACTOR}
git config user.email ${GITHUB_ACTOR}@users.noreply.github.com
git checkout gh-pages
git merge edit
## run the test script
cd /home
PYTHONPATH=/home python3 Psience/ci/tests/run_tests.py -v -d

# build docs and push
PYTHONPATH=/home python3 Psience/ci/build_docs.py
cd Psience
git add -A && git commit -m "Built out docs"
git push -u "https://$GITHUB_ACTOR:$GITHUB_TOKEN@github.com/McCoyGroup/Psience.git" gh-pages