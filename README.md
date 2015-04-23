# Network Evolution
Network Evolution and DMIs. 


To build: 
#g++ -o simulation -std=c++11 # -O2 -larmadillo main.cp
g++ -o simulation -std=c++11 -O2 -I -DARMA_DONT_USE_WRAPPER -framework Accelerate main.cp

To simulate:
./simulation num_pops num_individuals initial_x initial_y initial_reg_pattern theta gamma model num_generations allelic_stdev

Example:
./simulation 3 100 200 300 1122 300 1 'B' 100 200


For random number generation C++11 is required


Note: The armadillo library is no longer needed. 
For multiplying matrices the Armadillo library is used. Download library at
http://arma.sourceforge.net/download.html and follow instructions there to install. 

Once installed. Place all the contents of the include directory in the same directory as main.cp
