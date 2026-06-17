#!/usr/bin/env bash

if [[ "$VIRTUAL_ENV" == "" ]]; then
  echo "ERROR: Create and activate a virtual environment and try again"
  return
fi

pip install -e .

# Add TIMBERPATH to activate script if not already there
activate_path=$VIRTUAL_ENV/bin/activate

# 1. Use single quotes around the search pattern so grep looks for the literal '$PWD'
if ! grep -q 'export TIMBERPATH=\$PWD' "$activate_path"; then
  # 2. Use single quotes here to append the literal string instead of expanding it
  echo 'export TIMBERPATH="$PWD/TIMBER/"' >> "$activate_path"
  echo "Added TIMBERPATH to $activate_path"
fi

# Clone and build libarchive if not present
if [ ! -d "bin/libarchive" ]; then
  git clone -b v3.6.2 https://github.com/libarchive/libarchive.git
  cd libarchive
  cmake . -DCMAKE_INSTALL_PREFIX=../bin/libarchive
  make
  make install
  cd ..
  rm -rf libarchive
fi

if [ ! -f bin/restframes/lib/librestframes.so ]
then
  echo "Building RestFrames"
  cd bin/restframes/
  make
  cd ../../
fi

# Build libtimber if not present
if [ ! -d "bin/libtimber" ]; then
  make
fi
