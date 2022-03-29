#!/bin/bash

# Build iSCAM from source. If -c argument is present, call clean_admb script first.

# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
POSITIONAL=()
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -c|--clean)
      CLEAN=1
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL[@]}" # restore positional parameters

cwd=$(pwd)

cd /usr/bin/gfiscam
if [[ -n "${CLEAN}" ]]; then
  make clean
fi
make -j 7 dist

cd $cwd
