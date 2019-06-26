########################################################################
#
#
########################################################################

PROGRAM = polymer

OBJS = \
 main.o\
 inout.o\
 polymer.o\
 params.o\
 kdtree.o\
 montecarlo.o\
 energy.o

# The compiler used
CC = gcc

# List of external libraries used
LIBS += -lhdf5
LIBS += -lm

# The directory(s) where external libraries are storred
LDFLAGS += -L/usr/lib64
# The directory(s) where header files for external libraries are storred
INC += -I/usr/include 

# General compilation rules
CFLAGS += -Wall
CFLAGS += -O0
CFLAGS += -ggdb
#CFLAGS += -O3

# Rule to link the program
$(PROGRAM): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(LIBS) $^ -o $@

# Rule to compile the object files
%.o: %.c 
	$(CC) -c $(CFLAGS) $(INC) $< -o $@

.PHONY: clean

# Rule to clean all compiled objexts
clean:
	rm -f $(OBJS)
	rm -f $(PROGRAM)
