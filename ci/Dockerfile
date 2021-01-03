# GitHub Actions container to support doing CI
FROM python:3.8-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends git \
    && apt-get purge -y --auto-remove \
    && rm -rf /var/lib/apt/lists/*

# installs python packages we'd want for unit testing
RUN pip install \
              numpy \
              matplotlib \
              scipy \
              h5py

LABEL "name"="McCoyGroup/PsienceCI"
LABEL "maintainer"="Mark Boyer <b3m2a1uw.edu>"
LABEL "version"="0.0.1"

LABEL "com.github.actions.name"="Psience Continuous Integration"
LABEL "com.github.actions.description"="Runs CI for the Psience package"
LABEL "com.github.actions.icon"="git"
LABEL "com.github.actions.color"="orange"

RUN \
    cd /home && \
    git clone https://github.com/McCoyGroup/Peeves.git && \
    git clone -b edit https://github.com/McCoyGroup/McUtils.git && \
    git clone https://github.com/McCoyGroup/Psience.git

RUN \
    cd /home/Psience && \
    git checkout edit



ENTRYPOINT ["/bin/bash", "/home/Psience/ci/run_ci.sh"]