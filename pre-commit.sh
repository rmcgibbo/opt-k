#!/bin/sh
echo "Rebuilding README.md from LaTeX"
(cd "$(git rev-parse --show-toplevel)/paper" && make)
#cd "$(git rev-parse --show-toplevel)" || exit
touch .commit