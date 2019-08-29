.KEEP_STAT:

all: main

Compiler	= g++
FLAGS		= -D NDEBUG -O3 -m64
LIB		= -lm -lpthread -lstdc++
SOURCE		= main.cpp seq_clustering.cpp nr_db.cpp taxonomy.cpp tools.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

all:		main

main:		$(OBJECT)
			$(Compiler) $(FLAGS) $(OBJECT) -o MetaDB $(LIB)

%.o:		%.cpp $(HEADER)
			$(Compiler) $(FLAGS) -c $<
clean:
		rm -f *.o *~
