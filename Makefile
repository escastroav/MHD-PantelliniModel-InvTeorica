all: pantellini.x

%.o: %.cpp
	g++ -c $< -o $@

pantellini.x: pantellini.o main_pantellini.o
	g++ -std=c++11 -Wall -fconcepts -fsanitize=address -fsanitize=undefined $^ -o $@

test: pantellini.o model_tests.o
	g++ -std=c++11 -Wall -fsanitize=address -fsanitize=undefined $^ -o model_tests.x

valgrind: pantellini.x
	valgrind --tool=memcheck --leak-check=yes ./pantellini.x

clean:
	-rm -f *~ *# *o *out *x 
