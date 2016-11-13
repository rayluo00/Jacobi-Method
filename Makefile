LIBS=-lpthread

make: jacobi.c
	$(CC) -g -c jacobi.c
	$(CC) -g -o jacobi jacobi.o $(LIBS)

clean:
	$(RM) jacobi , *# , *~ , *.o , *.gch
