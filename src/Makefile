
CC = g++

GUROBIROOTDIR = /opt/gurobi811/linux64

CPLEXROOTDIR = /opt/cplex_12.9.0/CPLEX_Studio

CPLEXDIR = ${CPLEXROOTDIR}/cplex

CONCERTDIR = ${CPLEXROOTDIR}/concert

SUFSORTDIR = ./parallel-divsufsort

RHHDIR = ./robin-hood-hashing/src/include

CPPFLAGS = -I$(SUFSORTDIR)/include -I$(RHHDIR)

CPPFLAGS += -L$(SUFSORTDIR)/lib -L$(SUFSORTDIR)/external/libprange/lib -Wall -O2 -std=c++11

CPPFLAGS_CPLEX = -DCPLEX -D_GLIBCXX_DEBUG -DIL_STD

CPPFLAGS_CPLEX += -I$(CONCERTDIR)/include -I$(CPLEXDIR)/include

CPPFLAGS_CPLEX += -L$(CONCERTDIR)/lib/x86-64_linux/static_pic -L$(CPLEXDIR)/lib/x86-64_linux/static_pic

CPPFLAGS_GUROBI = -DGUROBI -I$(GUROBIROOTDIR)/include -L$(GUROBIROOTDIR)/lib
 
LIBS = -lpthread -fopenmp -ldivsufsort -llibprange

LIBS_CPLEX = -lconcert -lilocplex -lcplex -ldl

LIBS_GUROBI = -lgurobi_g++5.2 -lgurobi81

SRCS = gsa.cpp hashtrie.cpp binaryio.cpp build.cpp query.cpp main.cpp

EXEC = ../cammiq

cplex:
	$(CC) -o $(EXEC) $(SRCS) $(CPPFLAGS) $(CPPFLAGS_CPLEX) $(LIBS) $(LIBS_CPLEX)

gurobi:
	$(CC) -o $(EXEC) $(SRCS) $(CPPFLAGS) $(CPPFLAGS_GUROBI) $(LIBS) $(LIBS_GUROBI)
 
clean: 
	rm -f ./*.o $(EXEC) ../*.out

downloads:
	git clone https://github.com/jlabeit/parallel-divsufsort
	
downloadr:
	git clone https://github.com/martinus/robin-hood-hashing

divsufsort:
	cd $(SUFSORTDIR) && cmake . -DOPENMP=ON -DCILKP=OFF && make

