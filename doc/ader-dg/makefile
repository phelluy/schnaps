CC=gcc
OPTIMS=-std=c99 -lm  -O
LIBS= 
go: aderdg.o 
	$(CC) aderdg.o -o go $(OPTIMS) $(LIBS)

aderdg.o: aderdg.c aderdg.h
	$(CC) -c aderdg.c $(OPTIMS)

.PHONY clean: 
	 rm *.o go 
