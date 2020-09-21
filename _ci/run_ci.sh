#!/bin/bash
set -e

##
## Runs CI
##

# get into the parent folder
cd /home
# configure git in case we want to push stuff in the future
git config --global user.name ${GITHUB_ACTOR}
git config --global user.email ${GITHUB_ACTOR}@users.noreply.github.com
# clone in Peeves and McUtils for Unit testing
git clone https://github.com/McCoyGroup/Peeves.git
git clone https://github.com/McCoyGroup/McUtils.git
## run the testing script
cd /home
PYTHONPATH=/home python3 Psience/_ci/tests/run_tests.py -v -d
## write some artifacts, maybe, in the future (e.g. generated data from data_gen tests or images from matplotlib)