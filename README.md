# Network Evolution
Network Evolution and DMIs. 


To build: 
g++ -o simulation -std=c++11 main.cp

To simulate:
./simulation -p <num_pops> -i <num_individuals> -x <initial_x> -y <initial_y> -n <initial_network_pattern> -t <theta>  -g <gamma> -m <model_A_or_B> -G <num_generations> -d <allelic_stdev> -X <x_opt> -Y <y_opt> -s <omega11> -S <omega12> -r <percent_recombination>

Example:
./simulation -p 2 -i 5 -x 300 -y 300 -n 1122 -t 300 -g 1 -m B -G 20 -d 1 -X 300 -Y 300 -s 10000 -S 5000 -r 50

Notes:
For random number generation C++11 is required


Developmental tools for troubleshooting memory leaks:
valgrind --leak-check=full --show-leak-kinds=all ./simulation -p 2 -i 5 -x 300 -y 300 -n 1122 -t 300 -g 1 -m B -G 20 -d 1 -X 300 -Y 300 -s 10000 -S 5000 -r 50


