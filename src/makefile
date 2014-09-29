CC=gcc
OPTIMS=-std=c99 
LIBS=-lm 
go: geometry.o macromesh.o main.o 
	$(CC) geometry.o macromesh.o main.o -o go $(OPTIMS) $(LIBS)

geometry.o: geometry.c geometry.h global.h h20phi.h
	$(CC) -c geometry.c $(OPTIMS)

macromesh.o: macromesh.c macromesh.h global.h
	$(CC) -c macromesh.c $(OPTIMS)

	$(CC) -c main.c $(OPTIMS)

.PHONY clean: 
	 rm *.o ;  rm go 