
CC = g++ 

CPPFLAGS = -I./parallel-divsufsort/include -L./parallel-divsufsort/lib -L./parallel-divsufsort/external/libprange/lib -Wall -O2 -std=c++11
 
LIBS = -lpthread -fopenmp -ldivsufsort -llibprange

SRCS = gsa.cpp hashtrie.cpp binaryio.cpp build.cpp main.cpp
 
OBJS = $(SRCS: .cpp = .o) 

EXEC = build_index

$(EXEC) : $(OBJS) 
	$(CC) -o $@ $^ $(CPPFLAGS) $(LIBS)

	#$(CC) $(CPPFLAGS) $(LIBS) $^ -o $@  
 
%.o : %.cpp 
	$(CC) -c $(CPPFLAGS) $<  
 
clean: 
	rm -f *.o $(EXEC) *.out 
