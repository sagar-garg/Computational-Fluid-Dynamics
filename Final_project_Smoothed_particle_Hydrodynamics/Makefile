CC = gcc
CFLAGS = -Wall -pedantic 
.c.o:  ; $(CC) -c $(CFLAGS) $<

OBJ = 	helper.o\
      	main.o\
		init.o\
		uvp.o\
		visual.o\
      
all:  $(OBJ)
	$(CC) $(CFLAGS) -o sim $(OBJ)  -lm

%.o : %.c
	$(CC) -c $(CFLAGS) $*.c -o $*.o

clean:
	rm $(OBJ) sim

helper.o      : helper.h 
init.o        : helper.h init.h
uvp.o         : helper.h init.h uvp.h
visual.o      : visual.h helper.h


main.o        : helper.h init.h uvp.h visual.h

