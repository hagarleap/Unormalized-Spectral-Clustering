CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors -lm

# Specify the target executable and the source files needed to build it
my_app: spkmeans.o spkmeans.h
	$(CC) -o my_app spkmeans.o $(CFLAGS) 

# Specify the object files that are generated from the corresponding source files
spkmeans.o: spkmeans.c
	$(CC) -c spkmeans.c $(CFLAGS)



