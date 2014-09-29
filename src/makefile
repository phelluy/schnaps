CC=gcc
OPTIMS=-std=c99 
LIBS=-lm 
go: h20.o macromesh.o main.o 
	$(CC) h20.o macromesh.o main.o -o go $(OPTIMS) $(LIBS)

h20.o: h20.c h20.h global.h h20phi.h
	$(CC) -c h20.c $(OPTIMS)

macromesh.o: macromesh.c macromesh.h global.h
	$(CC) -c macromesh.c $(OPTIMS)

main.o: main.c macromesh.h global.h
	$(CC) -c main.c $(OPTIMS)

.PHONY clean: 
	 rm *.o ;  rm go 