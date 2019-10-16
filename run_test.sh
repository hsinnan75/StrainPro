#!/bin/bash
echo
echo "Test1 -- Download taxonomy files"
echo "Command=./download_taxonomy.sh"
echo
./download_taxonomy.sh

echo
echo "Test2 -- Build a reference database"
echo "Command=bin/StrainPro-build -r test/ecoli.fa -o ecoli_db"
echo
bin/StrainPro-build -r test/ecoli.fa -o ecoli_db

echo
echo "Test3 -- Read mapping & taxonomic composition analysis"
echo "Command=bin/StrainPro-map -i ecoli_db/ -f test/read.fq.gz"
echo
bin/StrainPro-map -i ecoli_db/ -f test/read.fq.gz

echo
rm -rf ecoli_db
echo "[End of test]"
