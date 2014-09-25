CC=gcc
OPTIMS= -g
LIBS= -lm  
go: macromesh.o main.o 
	$(CC) macromesh.o main.o -o go $(OPTIMS) $(LIBS)

macromesh.o: macromesh.c macromesh.h global.h
	$(CC) -c macromesh.c $(OPTIMS)

main.o: main.c macromesh.h global.h
	$(CC) -c main.c $(OPTIMS)

.PHONY clean: 
	 rm *.o rm go 
