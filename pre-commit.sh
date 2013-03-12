#!/bin/sh
files=`git diff --cached --name-status`
re="<files of importance>"
if [[ $files =~ $re ]]
then
  echo "Creating files"
  bundle exec create_my_files
  git add my_files
fi

exit 1