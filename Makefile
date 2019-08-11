all: get_cuts

get_cuts: get_cuts.cpp LineReader.h Progress.h
	clang++ -std=c++17 -Wall -O3 -g get_cuts.cpp -o $@
