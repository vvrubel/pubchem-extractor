#!/bin/sh

# Run script from project directory that contains `.env` file

file=$JSON
dir=${OUT_DIR}
f_dir=${dir}/files
name=$MONGO_DB_COLLECTION

rm -rf $dir

echo Downloading data into $dir
poetry run python -m molecad.data.console fetch \
       --out-dir $dir \
       --start 1 \
       --stop 1001 \
       --size 100
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
