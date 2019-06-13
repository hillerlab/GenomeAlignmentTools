KENTSRC_LIBS = $(KENTSRC_DIR)/lib/$(MACHTYPE)/jkweb.a $(KENTSRC_DIR)/htslib/libhts.a 

INCS = -I$(KENTSRC_DIR)/inc
LIBS = $(KENTSRC_LIBS) -lz -lssl -lcrypto -pthread

CC=gcc
CFLAGS= ${INCS} -O4 -static -Wall -Wformat -Wimplicit -Wreturn-type -Wuninitialized
LDFLAGS=

OBJECTS = scoreChain.o 

BIN=../../bin

scoreChain: $(OBJECTS)
	$(CC) $(OBJECTS) -o $(BIN)/$@ $(LDFLAGS) $(LIBS)

all: scoreChain

clean:
	rm -f *.o scoreChain
