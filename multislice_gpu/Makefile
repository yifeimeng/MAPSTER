# Edit the following variables as needed
INSTALL_PATH = /usr/local/cuda
EXTLIB = -L$(INSTALL_PATH)/lib64
CC	= nvcc
CFLAGS	= -x cu -Wno-deprecated-gpu-targets

INCLUDE	= -I$(INSTALL_PATH)/include
LIBSHDF = $(EXTLIB)

test: mslib.o misc.o main.o  
	$(CC) $? -o $@ $(LIBSHDF) -lcuda -lcudart 
main.o:
	$(CC) $(CFLAGS) -c main.c -o $@ $(INCLUDE)

mslib.o: 
	$(CC) $(CFLAGS) -c mslib.c -o $@ $(INCLUDE)

misc.o:
	$(CC) $(CFLAGS) -c misc.c -o $@ $(INCLUDE)
	 

clean: 
	rm *.o test
