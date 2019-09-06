get_cuts: *.cpp *.h
	clang++ -std=c++17 -Wall -O3 -g *.cpp -o $@
	./get_cuts --test
