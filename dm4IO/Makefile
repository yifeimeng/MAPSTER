# Edit the following variables as needed
INSTALL_PATH = /usr/local/hdf5
EXTLIB = -L$(INSTALL_PATH)/lib
CC	= gcc
CFLAGS	= 

INCLUDE	= -I$(INSTALL_PATH)/include
LIBSHDF = $(EXTLIB)

dm2h5: dm2h5.o parseDMx.o
	$(CC) $? -o $@ $(LIBSHDF) -lhdf5

sample: parseDMx.o sample.o 
	$(CC) $? -o $@

sample.o:
	$(CC) $(CFLAGS) -c sample.c -o $@ $(INCLUDE)

dm2h5.o: 
	$(CC) $(CFLAGS) -c dm2h5.c -o $@ $(INCLUDE)

parseDMx.o:
	$(CC) $(CFLAGS) -c parseDMx.c -o $@ $(INCLUDE)

	 

clean: 
	rm *.o test
