#!/bin/bash


echo Downloading data into ./data/fetch
poetry run python -m molecad.data.console fetch \
       --start 700001 \
       --stop 800001 \
       --req-size 100 \
       --f-size 1000


echo Split $file to $f_dir
poetry run python -m molecad.data.console split \
       --file $file \
       --f-dir $f_dir \
       --size 1000

echo Import files from directory $f_dir to MongoDB
poetry run python -m molecad.data.console populate \
       --f-dir $f_dir \
       --collection $name
