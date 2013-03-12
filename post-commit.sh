#!/bin/sh
#http://stackoverflow.com/questions/3284292/can-a-git-hook-automatically-add-files-to-the-commit
set -e
if [ -a "$(git rev-parse --show-toplevel)/.commit" ]
then
    rm "$(git rev-parse --show-toplevel)/.commit"
    git add paper/README.md
    git commit --amend -C HEAD --no-verify
    echo "Adding new README.md"
fi