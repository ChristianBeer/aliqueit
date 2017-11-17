#
# aliqueit Makefile for linux
#

# you'll need to change this if your gmp headers aren't installed globally:
GMP_DIR = ../../../../gmp

CC = g++
CFLAGS = -O3 -I$(GMP_DIR)

# you'll need to change this if your gmp headers aren't installed globally:
LIBS = -L/usr/local/lib -lgmp

OBJ = aliqueit.o cfg.o misc.o


aliqueit: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LIBS) -o aliqueit

%.o: %.cc
	$(CC) $(CFLAGS) -c -o $@ $<
clean:
	rm -f $(OBJ)
