CC=gcc
OPTIMS=-std=c99 
LIBS=-lm 
go: geometry.o interpolation.o macromesh.o main.o model.o test.o 
	$(CC) geometry.o interpolation.o macromesh.o main.o model.o test.o -o go $(OPTIMS) $(LIBS)

geometry.o: geometry.c geometry.h global.h h20phi.h
	$(CC) -c geometry.c $(OPTIMS)

interpolation.o: interpolation.c interpolation.h global.h
	$(CC) -c interpolation.c $(OPTIMS)

macromesh.o: macromesh.c macromesh.h global.h
	$(CC) -c macromesh.c $(OPTIMS)

main.o: main.c macromesh.h global.h geometry.h test.h
	$(CC) -c main.c $(OPTIMS)

model.o: model.c model.h
	$(CC) -c model.c $(OPTIMS)

test.o: test.c macromesh.h global.h geometry.h interpolation.h test.h \
  model.h
	$(CC) -c test.c $(OPTIMS)

.PHONY clean: 
	 rm *.o ;  rm go 