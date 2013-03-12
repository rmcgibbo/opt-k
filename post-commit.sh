#!/bin/sh
#http://stackoverflow.com/questions/3284292/can-a-git-hook-automatically-add-files-to-the-commit
if [ -a .commit ]
    then
    cd "$(git rev-parse --show-toplevel)/paper"
    del .commit
    git add README.md
    git commit --ammend -C HEAD --no-verify
fi
exit