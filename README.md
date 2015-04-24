# Network Evolution
Network Evolution and DMIs. 


To build: 
g++ -o simulation -std=c++11 main.cp

To simulate:
./simulation num_pops num_individuals initial_x initial_y initial_reg_pattern theta gamma model num_generations allelic_stdev x_opt y_opt omega11 omega12

Example:
./simulation 3 100 200 300 1122 300 1 'B' 100 200 300 300 10000 0

Notes:
For random number generation C++11 is required


Developmental tools for troubleshooting memory leaks:
valgrind --leak-check=full --show-leak-kinds=all ./simulation 30 1000 200 300 1122 300 1 'B' 100 200 300 300 10000 0


