# Edit the following variables as needed
INSTALL_PATH =
EXTLIB = 
CC	= gcc
CFLAGS	= 

INCLUDE	= 
LIBSHDF = 

test: crystal.o main.o  
	$(CC) $? -o $@

main.o: main.c
	$(CC) $(CFLAGS) -c $? -o $@ 

crystal.o: crystal.c
	$(CC) $(CFLAGS) -c $? -o $@ 

	 
clean: 
	rm *.o test
