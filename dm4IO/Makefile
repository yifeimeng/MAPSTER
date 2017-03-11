# Edit the following variables as needed
INSTALL_PATH = /usr/local/hdf5
EXTLIB = -L$(INSTALL_PATH)/lib
CC	= gcc
CFLAGS	= 

INCLUDE	= -I$(INSTALL_PATH)/include
LIBSHDF = $(EXTLIB)

test: parseDMx.o main.o 
	$(CC) $? -o $@ $(LIBSHDF) -lhdf5 
main.o:
	$(CC) $(CFLAGS) -c main.c -o $@ $(INCLUDE)

pinpointTagsDM4.o:
	$(CC) $(CFLAGS) -c parseDMx.c -o $@ $(INCLUDE)
	 

clean: 
	rm *.o test
