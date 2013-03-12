#!/bin/sh
echo "Rebuilding README.md from LaTeX"
cd "$(git rev-parse --show-toplevel)/paper"
make
git add README.md
touch .commit