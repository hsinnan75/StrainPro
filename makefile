.KEEP_STAT:

all:		index nr db map

index:
		make -C src/BWT_Index && mv src/BWT_Index/bwt_index .

nr:
		make -C src/MetaNR && mv src/MetaNR/MetaNR StrainPro-nr
		
db:
		make -C src/MetaDB && mv src/MetaDB/MetaDB StrainPro-build

map:
		make -C src/MetaMapping && mv src/MetaMapping/MetaMapping StrainPro-map

clean:
		make clean -C src/BWT_Index
		make clean -C src/MetaNR
		make clean -C src/MetaDB
		make clean -C src/MetaMapping
