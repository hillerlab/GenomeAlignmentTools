KENTSRC_LIBS = $(KENTSRC_DIR)/lib/$(MACHTYPE)/jkhgap.a  $(KENTSRC_DIR)/lib/$(MACHTYPE)/jkweb.a $(KENTSRC_DIR)/htslib/libhts.a 

INCS = -I$(KENTSRC_DIR)/inc -I$(KENTSRC_DIR)/hg/inc
LIBS = $(KENTSRC_LIBS) -lm  -lz -lcrypt -lnsl -lssl -lcrypto -lstdc++ -lrt -pthread

CC=gcc
CFLAGS= ${INCS} -O4 -static -Wall -Wformat -Wimplicit -Wreturn-type -Wuninitialized
LDFLAGS=

OBJECTS = chainCleaner.o 

BIN=../../bin

chainCleaner: $(OBJECTS)
	$(CC) $(OBJECTS) -o $(BIN)/$@ $(LDFLAGS) $(LIBS)

all: chainCleaner

clean:
	rm -f *.o chainCleaner
