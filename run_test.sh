#!/bin/bash
echo
echo "Test1 -- Download taxonomy files"
echo "Command=./download_taxonomy.sh"
echo
./download_taxonomy.sh

echo
echo "Test2 -- Build a reference database"
echo "Command=./StrainPro-build -r test/ecoli.fa -o ecoli_db"
echo
./StrainPro-build -r test/ecoli.fa -o ecoli_db

echo
echo "Test3 -- Read mapping & taxonomic composition analysis"
echo "Command=./StrainPro-map -i ecoli_db/ -f test/read.fq.gz"
echo
./StrainPro-map -i ecoli_db/ -f test/read.fq.gz

echo
echo "[End of test]"
