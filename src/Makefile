#
# aliqueit Makefile for linux
#

# you'll need to change this if your gmp headers are installed globally:
# you can also place symlinks in ../3rdParty/ to point to your prefered GMP version
GMP_INCLUDE_DIR = ../3rdParty/include
GMP_LIBRARY_DIR = ../3rdParty/lib

CC = g++
CFLAGS = -O3 -I$(GMP_INCLUDE_DIR)

LIBS = -L$(GMP_LIBRARY_DIR) -lgmp

OBJ = aliqueit.o cfg.o misc.o


aliqueit: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LIBS) -o aliqueit

%.o: %.cc
	$(CC) $(CFLAGS) -c -o $@ $<
clean:
	rm -f $(OBJ)
