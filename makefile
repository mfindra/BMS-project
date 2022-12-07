AUTHOR= xfindr00
NAME= bms
CC= g++
FLAGS= -std=c++11 -pedantic

all:
	$(CC) $(NAME).cpp -o $(NAME) $(FLAGS)  

clean:
	rm $(NAME)

pack: 
	tar -cf $(AUTHOR).tar manual.pdf bms.cpp makefile bms.1