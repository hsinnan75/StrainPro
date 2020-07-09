.KEEP_STAT:

all: main

CXX		= g++
FLAGS		= -Wall -D NDEBUG -O3 -m64
LIB		= -lm -lz -lpthread -lstdc++
SOURCE		= main.cpp GetData.cpp Mapping.cpp Typing.cpp ReadAlignment.cpp taxonomy.cpp nw_alignment.cpp tools.cpp bwt_index.cpp bwt_search.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

all:		main

main:		$(OBJECT)
			$(CXX) $(FLAGS) $(OBJECT) -o StrainPro-map $(LIB)

%.o:		%.cpp $(HEADER)
			$(CXX) $(FLAGS) -c $<
clean:
		rm -f *.o *~
