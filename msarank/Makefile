# Makefile for msarank

LIB = -lm -lz
CC = gcc
CFLAGS = -O -Wall -Werror

# compiler flags and timings
# opt   user  sys
# none  58.02 43.84
# -O    18.33 10.53
# -O1   14.48 14.82
# -O2   21.06 23.60
# -O3   17.03 20.29

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


