#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

SOURCE_BRANCH="master"
TARGET_BRANCH="gh-pages"

function doCompile {
  cp -ra dox/html/* out
}

# Only first worker builds documentation
# echo "TRAVIS_JOB_NUMBER=$TRAVIS_JOB_NUMBER"
# echo "WORKER=${TRAVIS_JOB_NUMBER: -1}"
# if [ "${TRAVIS_JOB_NUMBER: -1}" != "1" ]; then
#     echo "Skipping deploy; Not my job."
#     exit 0
# fi


# Pull requests and commits to other branches shouldn't try to deploy, just build to verify
if [ "$TRAVIS_PULL_REQUEST" != "false" ]; then
    echo "Skipping deploy; This is a pull request."
    exit 0
fi
if [ "$TRAVIS_BRANCH" != "$SOURCE_BRANCH" ]; then
    echo "Skipping deploy; This is not the main branch."
    exit 0
fi

# Save some useful information
REPO=`git config remote.origin.url`
SSH_REPO=${REPO/https:\/\/github.com\//git@github.com:}
SHA=`git rev-parse --verify HEAD`

# Clone the existing gh-pages for this repo into out/
# Create a new empty branch if gh-pages doesn't exist yet (should only happen on first deply)
git clone $REPO out
cd out
git checkout $TARGET_BRANCH || git checkout --orphan $TARGET_BRANCH
cd ..

# Clean out existing contents
rm -rf out/* || exit 0

# Run our compile script
doCompile

# Now let's go have some fun with the cloned repo
cd out
git config user.name "Travis CI"
git config user.email "$COMMIT_AUTHOR_EMAIL"

# If there are no changes to the compiled out (e.g. this is a README update) then just bail.
if [ -z `git diff --exit-code` ]; then
    echo "No changes to the output on this push; exiting."
    exit 0
fi

# Commit the "changes", i.e. the new version.
# The delta will show diffs between new and old versions.
touch .nojekyll
git add .nojekyll
## _main.pdf keeps changing - probably time related. For now, get rid of it:
rm -f _book/_main.pdf
git add --all .
git commit -m "Deploy to GitHub Pages: ${SHA}"

# Get the deploy key by using Travis's stored variables to decrypt deploy_key.enc
openssl aes-256-cbc -K $encrypted_6397b83665d1_key -iv $encrypted_6397b83665d1_iv -in ../bin/deploy_key.enc -out ../bin/deploy_key -d
chmod 600 ../bin/deploy_key
eval `ssh-agent -s`
ssh-add ../bin/deploy_key

# Now that we're all set up, we can push.
git push $SSH_REPO $TARGET_BRANCH
