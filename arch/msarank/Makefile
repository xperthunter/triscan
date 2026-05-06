# Makefile for msarank

LIB = -lm -lz
CC = gcc
CFLAGS = -O3 -Wall -Werror

# compiler flags and timings
#       Linux-x86-VM  MacOS-ARM64
# opt   user  sys     user  sys
# none  58.02 43.84   92.22 91.56
# -O    18.33 10.53   22.89 22.55
# -O1   14.48 14.82   22.88 22.41
# -O2   21.06 23.60   17.22 16.88
# -O3   17.03 20.29   14.60 14.15

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


