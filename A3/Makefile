CC=mpic++
CFLAGS=-std=c++17 -fopenmp -O2 -g

a3: main.o
	$(CC) $(CFLAGS) main.o -o a3

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp -o main.o

clean:
	rm a2 main.o
# mpirun -np 2 ./a3 --taskid=1 --inputpath=/test4/test-input-4.gra --headerpath=/test4/test-header-4.dat --outputpath=output.txt --verbose=1 --startk=1 --endk=8