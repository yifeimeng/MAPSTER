# Edit the following variables as needed
INSTALL_PATH = /usr/local/hdf5
EXTLIB = -L$(INSTALL_PATH)/lib
CC	= gcc
CFLAGS	= 

INCLUDE	= -I$(INSTALL_PATH)/include
LIBSHDF = $(EXTLIB)

test: pinpointTagsDM4.o main.o 
	$(CC) $? -o $@ 
main.o:
	$(CC) $(CFLAGS) -c main.c -o $@ $(INCLUDE)

pinpointTagsDM4.o:
	$(CC) $(CFLAGS) -c pinpointTagsDM4.c -o $@ $(INCLUDE)
	 

clean: 
	rm *.o test
