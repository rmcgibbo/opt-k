#!/bin/sh
set -e
echo "Rebuilding README.md from LaTeX"
(cd "$(git rev-parse --show-toplevel)/paper" && make)
touch .commit