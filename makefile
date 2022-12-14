#AUTHOR= xfindr00
#NAME= bms
#CC= g++
#FLAGS= -std=c++11 -pedantic

#all:
#	$(CC) $(NAME).cpp -o $(NAME) $(FLAGS)  

#clean:
#	rm $(NAME)

#pack: 
#	tar -cf $(AUTHOR).tar manual.pdf bms.cpp makefile bms.1


# autgor for packing
AUTHOR = 221646

# program name
NAME = bms

# define the C++ compiler to use
CC = g++

# define any compile-time flags
CFLAGS = -std=c++11 -pedantic

# define the C++ source files
SRCS = main.cc rs.cc gf.cc

# define the C++ object files
OBJS = $(SRCS:.cc=.o)

# define the executable file
MAIN = main

#
# The default target, which compiles the executable
#

all: $(MAIN)

#
# The rule to link the object files to create the executable
#

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) -o $(NAME) $(OBJS)

#
# The rule to compile the C++ source files
#

.cc.o:
	$(CC) $(CFLAGS) -c $<

#
# A rule to clean up the directory
#

clean:
	rm -f $(MAIN) $(OBJS) $(NAME)
