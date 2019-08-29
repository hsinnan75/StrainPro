.KEEP_STAT:

all: main

Compiler	= g++
FLAGS		= -D NDEBUG -O3 -m64
LIB		= -lm -lpthread -lstdc++
SOURCE		= main.cpp SeqTaxing.cpp tools.cpp bwt_index.cpp bwt_search.cpp taxonomy.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

all:		main

main:		$(OBJECT)
			$(Compiler) $(FLAGS) $(OBJECT) -o StrainPro-nr $(LIB)

%.o:		%.cpp $(HEADER)
			$(Compiler) $(FLAGS) -c $<
clean:
		rm -f *.o *~
