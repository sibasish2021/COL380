all: exec

exec: main.o
	g++ -fopenmp -Wl,-rpath=. -g -o exec main.o -L. -l:library.so

main.o: main.cpp library.hpp
	g++ -fopenmp -c main.cpp

clean:
	rm -f *.o exec 

