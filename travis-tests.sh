#!/bin/bash

PATH=/opt/cms/common:/opt/cms/bin:$PATH
cmssw_version=${TRAVIS_BRANCH#forest_}

source /opt/cms/cmsset_default.sh

scramv1 project CMSSW $cmssw_version
cd $cmssw_version/src
eval `scramv1 runtime -sh`

git config --global user.name docker-tmp
git config --global user.email "cmsbuild@docker-tmp"

printf '%s\n' N | git cms-init --upstream-only
git cms-merge-topic -u bi-ran:$TRAVIS_BRANCH

git checkout -b travis-$TRAVIS_PULL_REQUEST_BRANCH
git pull --no-edit https://github.com/$TRAVIS_PULL_REQUEST_SLUG.git \
    $TRAVIS_PULL_REQUEST_BRANCH

pushd HeavyIonsAnalysis/JetAnalysis/python/jets
./makeJetSequences.sh
popd

scram b -j$(nproc)

pushd HeavyIonsAnalysis/JetAnalysis/test
./tests.sh -t -j $(nproc)
