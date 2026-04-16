# Makefile for het-vs-hom

LIB = -lm -lz
CC = gcc
CFLAGS = -O2 -Wall -Werror

APP = msarank
SRC = msarank.c
OBJ = msarank.o

# Targets

default:
	make $(APP)

$(APP): $(OBJ)
	$(CC) -o $(APP) $(OBJ) $(LIB)

clean:
	rm -f *.o $(APP)

# Inference Rules

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

