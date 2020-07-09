.KEEP_STAT:

all: main

CXX		= g++
FLAGS		= -Wall -D NDEBUG -O3 -m64
LIB		= -lm -lpthread -lstdc++
SOURCE		= main.cpp SeqTaxing.cpp tools.cpp bwt_index.cpp bwt_search.cpp taxonomy.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

all:		main

main:		$(OBJECT)
			$(CXX) $(FLAGS) $(OBJECT) -o StrainPro-rep $(LIB)

%.o:		%.cpp $(HEADER)
			$(CXX) $(FLAGS) -c $<
clean:
		rm -f *.o *~
