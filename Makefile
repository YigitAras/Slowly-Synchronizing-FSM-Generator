CC=g++
CFLAGS=-O3 -std=c++11
GEN_FILE=genetic_generator.cpp
GRE_FILE=greedy_generator.cpp

all: gen greedy

gen: $(GEN_FILE)
	$(CC) $(GEN_FILE) $(CFLAGS) -o Genetic_Gen
	@echo "Genetic is compiled" 

greedy: $(GRE_FILE)
	$(CC) $(GRE_FILE) $(CFLAGS) -o Greedy_Gen
	@echo "Greedy is compiled"
clean:
	rm -rf Greedy_Gen Genetic_Gen
