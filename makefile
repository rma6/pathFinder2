compile: main.cpp
	g++ -std=c++17 -m64 -g -pthread -fopenmp -O3 -Wall -Wextra main.cpp -o exe
