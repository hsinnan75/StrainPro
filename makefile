.KEEP_STAT:

all:		index rep db map

index:
		make -C src/BWT_Index && mv src/BWT_Index/bwt_index .

rep:
		make -C src/MetaNR && mv src/MetaNR/StrainPro-rep .
		
db:
		make -C src/MetaDB && mv src/MetaDB/StrainPro-build .

map:
		make -C src/MetaMapping && mv src/MetaMapping/StrainPro-map .

clean:
		make clean -C src/BWT_Index
		make clean -C src/MetaNR
		make clean -C src/MetaDB
		make clean -C src/MetaMapping
