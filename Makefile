#Compiler
CC = g++
FLAGS = -Wall -g 

#Source folder
SRC = $(realpath $(CURDIR))/src/

all: genotype_concordance

genotype_concordance: concordance.o command_line_parser.o
	$(CC) $(FLAGS) -o GenotypeConcordance  concordance.o command_line_parser.o

concordance.o: $(SRC)concordance.cpp
	$(CC) $(FLAGS) -c $(SRC)concordance.cpp

command_line_parser.o: $(SRC)command_line_parser.cpp $(SRC)command_line_parser.hpp
	$(CC) $(FLAGS) -c $(SRC)command_line_parser.cpp

clean:
	rm -f concordance.o command_line_parser.o GenotypeConcordance
