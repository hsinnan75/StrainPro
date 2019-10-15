.KEEP_STAT:

all: main

CXX		= g++
FLAGS		= -Wall -D NDEBUG -O3 -m64
LIB		= -lm -lpthread -lstdc++
SOURCE		= main.cpp seq_clustering.cpp nr_db.cpp taxonomy.cpp tools.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

all:		main

main:		$(OBJECT)
			$(CXX) $(FLAGS) $(OBJECT) -o StrainPro-build $(LIB)

%.o:		%.cpp $(HEADER)
			$(CXX) $(FLAGS) -c $<
clean:
		rm -f *.o *~
