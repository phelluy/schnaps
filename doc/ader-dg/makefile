CC=gcc
OPTIMS= -O
LIBS= -lm
go: *.o 
	$(CC) *.o -o go $(OPTIMS) $(LIBS)

aderdg.o: aderdg.c aderdg.h
	$(CC) -c aderdg.c $(OPTIMS)

.PHONY clean: 
	 rm *.o rm go 
