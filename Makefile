all: pantellini.x

%.o: %.cpp
	g++ -c $< -o $@

pantellini.x: pantellini.o main_pantellini.o
	g++ -fconcepts -std=c++11 -fsanitize=address -fsanitize=undefined $^ -o $@

test: main_project.o project_test.o
	g++ -std=c++11 -Wall -fsanitize=address -fsanitize=undefined $^ -o project_test.o -lgsl -lgslcblas

valgrind: main_project.x
	valgrind --tool=memcheck --leak-check=yes ./main_project.x

clean:
	-rm -f *~ *# *o *out *x 
