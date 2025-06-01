CFLAG = -std=c++17 -O2 -I./
LDFLAG = -lgflags_nothreads 
all: main #preprocess

main: main.cpp
	g++ -g -o main main.cpp graph.hpp $(CFLAG) $(LDFLAG)

preprocess: preprocess.cpp
	g++ -g -o preprocess preprocess.cpp $(CFLAG) $(LDFLAG)
