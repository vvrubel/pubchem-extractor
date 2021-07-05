#!/bin/sh

source ../../.env

rm -rf $OUT_DIR

echo Downloading data into $OUT_DIR
poetry run python -m molecad.data.console fetch \
       --out-dir $OUT_DIR \
       --start 1 \
       --stop 1001 \
       --size 100


echo Split $FILE to $F_DIR
poetry run python -m molecad.data.console split \
       --file $FILE \
       --f-dir $F_DIR \
       --size 1000

echo Import files from directory $F_DIR
poetry run python -m molecad.data.console populate\
       --f-dir $F_DIR