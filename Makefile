compiler = icpc
flags =
flags = -I. -I./eigen/ -lm	-fopenmp -g

headers = $(wildcard *.hpp)
sources = $(wildcard *.cpp)
objects = $(sources:.cpp=.o)

executables: solver

%.o: %.cpp $(headers)
	$(compiler) -c -o $@ $< $(flags)

solver: $(objects)
	$(compiler) -o $@ $^ $(flags)

clean:
	rm -f *.o	solver

run:
	OMP_NUM_THREADS=4 ./solver
