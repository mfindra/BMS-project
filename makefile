#author: Michal Findra, xfindr00
#project: Reed-Salomon encoder
#date: 12.12.2022


AUTHOR= 221646
NAME= bms
CC= g++
FLAGS= -std=c++11 -pedantic
SOURCES= main.cc reed_salomon.cc galois_field.cc

all:
	$(CC) $(SOURCES) -o $(NAME) $(FLAGS)  

clean:
	rm $(NAME)

pack:
	zip 221646.zip main.cc reed_salomon.cc reed_salomon.hh galois_field.hh galois_field.cc Makefile
