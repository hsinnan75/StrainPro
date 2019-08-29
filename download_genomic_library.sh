#!/bin/bash

# Copyright 2013-2018, Derrick Wood <dwood@cs.jhu.edu>
#
# Download specific genomic libraries.

set -u  # Protect against uninitialized vars.
set -e  # Stop on error

LIBRARY_DIR="database"
NCBI_SERVER="ftp.ncbi.nlm.nih.gov"
FTP_SERVER="ftp://$NCBI_SERVER"
RSYNC_SERVER="rsync://$NCBI_SERVER"
THIS_DIR=$PWD
library_name="$1"
ftp_subdir=$library_name
library_file="library.fna"

case $library_name in
  "archaea" | "bacteria" | "viral" | "fungi" | "human" | "protozoa")
    mkdir -p $LIBRARY_DIR/$library_name
    cd $LIBRARY_DIR/$library_name
    rm -f assembly_summary.txt
    remote_dir_name=$library_name
    if [ "$library_name" = "human" ]; then
      remote_dir_name="vertebrate_mammalian/Homo_sapiens"
    fi          
    if ! wget -q $FTP_SERVER/genomes/refseq/$remote_dir_name/assembly_summary.txt; then
      1>&2 echo "Error downloading assembly summary file for $library_name, exiting."
      exit 1
    fi
    if [ "$library_name" = "human" ]; then
      grep "Genome Reference Consortium" assembly_summary.txt > x
      mv x assembly_summary.txt
    fi
    rm -rf all/ library.f* manifest.txt rsync.err
    $THIS_DIR/rsync_from_ncbi.pl assembly_summary.txt
    ;;
  *)
    1>&2 echo "Unsupported library.  Valid options are: "
    1>&2 echo "  archaea bacteria viral fungi protozoa"
    exit 1
    ;;
esac
