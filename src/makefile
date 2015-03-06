UNAME = $(shell uname)
CC=gcc
CFLAGS=-Wall -g -std=c99
#Try to compile on Mac OS with Libxml2 installed through MacPorts
ifeq ($(UNAME), Darwin)
LFLAGS=-lxml2 -lz -liconv -lm
IFLAGS=-I /opt/local/include/libxml2
else
#Try to compile on Linux with Libxml2 installed through apt-get
ifeq ($(UNAME), Linux)
LFLAGS=-lxml2 -lz -lm
IFLAGS=-I /usr/include/libxml2
#Try to install on Windows with Libxml2 installed through Cygwin
else
LFLAGS=-lxml2 -lz -liconv -lm
IFLAGS=-I cygwin64\usr\include\libxml2
endif
endif

all: ebtel

ebtel: ebtel_functions_loop.o ebtel_functions_param.o ebtel_functions_util.o ebtel_functions_solvers.o ebtel_functions_heating.o ebtel_main.o 
	$(CC) $(CFLAGS) ebtel_functions_loop.o ebtel_functions_param.o ebtel_functions_util.o ebtel_functions_solvers.o ebtel_functions_heating.o ebtel_main.o -o ../bin/ebtel $(LFLAGS)

ebtel_functions_loop.o: ebtel_functions_loop.c
	$(CC) -c  $(CFLAGS) $(IFLAGS) ebtel_functions_loop.c 

ebtel_functions_param.o: ebtel_functions_param.c
	$(CC) -c  $(CFLAGS) $(IFLAGS) ebtel_functions_param.c

ebtel_functions_util.o: ebtel_functions_util.c
	$(CC) -c $(CFLAGS) $(IFLAGS) ebtel_functions_util.c

ebtel_functions_solvers.o: ebtel_functions_solvers.c
	$(CC) -c  $(CFLAGS) $(IFLAGS) ebtel_functions_solvers.c

ebtel_functions_heating.o: ebtel_functions_heating.c
	$(CC) -c  $(CFLAGS) $(IFLAGS) ebtel_functions_heating.c

ebtel_main.o: ebtel_main.c
	$(CC) -c $(CFLAGS) $(IFLAGS) ebtel_main.c

clean: 
	rm *.o ../bin/ebtel
