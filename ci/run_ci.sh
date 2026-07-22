#!/bin/sh
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
repo="https://$GITHUB_ACTOR:$GITHUB_TOKEN@github.com/McCoyGroup/Psience.git"
git checkout gh-pages
git pull
git merge $branch
git push -u $repo gh-pages
## run the test script
cd /home

for i in {1..3}; do
if [[ "$1" == "tests" ]]; then
  shift;
  run_tests=true
fi
if [[ "$1" == "docs" ]]; then
  shift;
  build_docs=true
fi
if [[ "$1" == "stubs" ]]; then
  shift;
  build_stubs=true
fi
done

if [[ "$run_test" == "true" ]]; then
  if [[ "$branch" == "master" ]]; then
    PYTHONPATH=/home python3 Psience/ci/tests/run_tests.py -v -d
  else
    PYTHONPATH=/home python3 Psience/ci/tests/run_tests.py -d
  fi
fi

if [[ "$build_docs" == "true" ]]; then
  if [[ "$branch" == "master" ]]; then
    # build docs and push
    PYTHONPATH=/home python3 Psience/ci/build_docs.py
    cp -r Psience/ci/docs Psience/
    rm -rf Psience/ci/docs/Psience
    rm Psience/ci/docs/_config.yml
    cd Psience
    git add -A
    git diff-index --quiet HEAD || git commit -m "Built out docs"
    git push -u $repo gh-pages
  fi
fi

if [[ "$build_stubs" == "true" ]]; then
  if [[ "$branch" == "master" ]]; then
    # build stubs and push
    PYTHONPATH=/home python3 Psience/ci/build_stubs.py
    cp -r Psience/ci/stubs Psience/
    rm -rf Psience/ci/stubs
#    rm McUtils/ci/docs/_config.yml
#    cd McUtils
#    git add -A
#    git diff-index --quiet HEAD || git commit -m "Built out docs"
#    git push -u $repo gh-pages
  fi
fi