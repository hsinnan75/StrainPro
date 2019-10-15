.KEEP_STAT:

all:		index rep db map

bin:
		mkdir -p bin/
		
index: bin
		$(MAKE) -C src/BWT_Index
		cp -f src/BWT_Index/bwt_index bin/

rep: bin
		$(MAKE) -C src/MetaNR
		cp -f src/MetaNR/StrainPro-rep bin/
		
db: bin
		$(MAKE) -C src/MetaDB
		cp -f src/MetaDB/StrainPro-build bin/

map: bin
		$(MAKE) -C src/MetaMapping
		cp -f src/MetaMapping/StrainPro-map bin/

clean:
		rm -f bin/bwt_index bin/StrainPro-rep bin/StrainPro-build bin/StrainPro-map
		$(MAKE) clean -C src/BWT_Index
		$(MAKE) clean -C src/MetaNR
		$(MAKE) clean -C src/MetaDB
		$(MAKE) clean -C src/MetaMapping
