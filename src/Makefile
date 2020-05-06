CC = g++ 

CPPFLAGS = -D_GLIBCXX_DEBUG -DIL_STD -I./parallel-divsufsort/include -I./cplex_12.9.0/CPLEX_Studio/concert/include -I./cplex_12.9.0/CPLEX_Studio/cplex/include -L./cplex_12.9.0/CPLEX_Studio/concert/lib/x86-64_linux/static_pic -L./parallel-divsufsort/lib -L./parallel-divsufsort/external/libprange/lib -L./cplex_12.9.0/CPLEX_Studio/cplex/lib/x86-64_linux/static_pic -Wall -O2 -std=c++11
 
LIBS = -lpthread -fopenmp -ldivsufsort -llibprange -lconcert -lilocplex -lcplex -ldl

SRCS = gsa.cpp hashtrie.cpp binaryio.cpp build.cpp query.cpp main.cpp
 
OBJS = $(SRCS: .cpp = .o) 

EXEC = run

$(EXEC) : $(OBJS) 
	$(CC) -o $@ $^ $(CPPFLAGS) $(LIBS)

	#$(CC) $(CPPFLAGS) $(LIBS) $^ -o $@  
 
%.o : %.cpp 
	$(CC) -c $(CPPFLAGS) $<  
 
clean: 
	rm -f *.o $(EXEC) *.out
