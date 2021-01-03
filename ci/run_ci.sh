#!/bin/bash
set -e

##
## Runs CI
##

branch=$(git rev-parse --abbrev-ref HEAD)
echo "Running tests on $branch"

# get into the parent folder and merge in changes from the master branch
cd /home/Psience
git config user.name ${GITHUB_ACTOR}
git config user.email ${GITHUB_ACTOR}@users.noreply.github.com
git checkout gh-pages
git pull
git merge $branch
## run the test script
cd /home

if [[ "$branch" = "master" ]]; then
  PYTHONPATH=/home python3 Psience/ci/tests/run_tests.py -v -d
else
  PYTHONPATH=/home python3 Psience/ci/tests/run_tests.py -d
fi

if [[ "$branch" = "edit" ]]; then
  # build docs and push
  PYTHONPATH=/home python3 Psience/ci/build_docs.py
  cd Psience
  git add -A && git commit -m "Built out docs"
  git push -u "https://$GITHUB_ACTOR:$GITHUB_TOKEN@github.com/McCoyGroup/Psience.git" gh-pages
fi