#!/bin/sh
#http://stackoverflow.com/questions/3284292/can-a-git-hook-automatically-add-files-to-the-commit
if [ -a "$(git rev-parse --show-toplevel)/.commit" ]
then
    rm "$(git rev-parse --show-toplevel)/.commit"
    git add paper/README.md
    git commit --amend -C HEAD --no-verify
    echo "post commit hook"
fi