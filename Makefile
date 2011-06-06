
LIBDIR=/usr/local/lib
LIBS=-lcunit


inverse: myinverse.c
	gcc -o myinverse myinverse.c -L$(LIBDIR) $(LIBS)

clean:
	rm myinverse