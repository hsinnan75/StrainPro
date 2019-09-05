#!/bin/bash

#set -u  # Protect against uninitialized vars.
#set -e  # Stop on error

TAX_DIR="taxonomy"
RSYNC_SERVER="rsync://ftp.ncbi.nlm.nih.gov"

mkdir -p "$TAX_DIR"
cd "$TAX_DIR"

function download_file() {
  file="$1"
  rsync --no-motd ${RSYNC_SERVER}${file} .
}

1>&2 echo -n "Downloading taxonomy data..."
download_file "/pub/taxonomy/taxdump.tar.gz"
1>&2 echo " done."

1>&2 echo -n "Decompressing taxonomy tree data..."
tar zxf taxdump.tar.gz
1>&2 echo " done."

1>&2 echo -n "Removing unnecessary files..."
rm -f taxdump.tar.gz citations.dmp delnodes.dmp division.dmp gc.prt gencode.dmp readme.txt
1>&2 echo " done."

