CC = g++ 

CPLEXROOTDIR = ./cplex_12.9.0/CPLEX_Studio

CPLEXDIR = ${CPLEXROOTDIR}/cplex

CONCERTDIR = ${CPLEXROOTDIR}/concert

SUFSORTDIR = ./parallel-divsufsort

RHHDIR = ./robin-hood-hashing/src/include

CPPFLAGS = -D_GLIBCXX_DEBUG -DIL_STD 

CPPFLAGS += -I$(SUFSORTDIR)/include -I$(CONCERTDIR)/include -I$(CPLEXDIR)/include -I$(RHHDIR)

CPPFLAGS += -L$(CONCERTDIR)/lib/x86-64_linux/static_pic -L$(CPLEXDIR)/lib/x86-64_linux/static_pic 

CPPFLAGS += -L$(SUFSORTDIR)/lib -L$(SUFSORTDIR)/external/libprange/lib -Wall -O2 -std=c++11
 
LIBS = -lpthread -fopenmp -ldivsufsort -llibprange -lconcert -lilocplex -lcplex -ldl

SRCS = gsa.cpp hashtrie.cpp binaryio.cpp build.cpp query.cpp main.cpp
 
OBJS = $(SRCS: .cpp = .o) 

EXEC = ../cammiq

$(EXEC) : $(OBJS)
	cd $(SUFSORTDIR)
	$(CC) -o $@ $^ $(CPPFLAGS) $(LIBS)

	#$(CC) $(CPPFLAGS) $(LIBS) $^ -o $@  
 
%.o : %.cpp 
	$(CC) -c $(CPPFLAGS) $<  
 
clean: 
	rm -f *.o $(EXEC) *.out

download:
	git clone https://github.com/jlabeit/parallel-divsufsort
	git clone https://github.com/martinus/robin-hood-hashing

divsufsort:
	cd $(SUFSORTDIR) && cmake . -DOPENMP=ON -DCILKP=OFF && make

