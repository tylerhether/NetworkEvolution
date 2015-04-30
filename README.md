# Network Evolution
Network Evolution and DMIs. 


I still need to add the following things:
* 1) migration
* 2) regulatory mutations
* 3) test for hybrid incompatibility 
* 4) Calculate other summary statistics and output them to file



To build: 
g++ -std=c++11 main.cp -o simulation

To simulate:
./simulation -npops <n_populations> -i <n_individuals> -gen <n_generations> -mu <coding_mutation_rate> -mu_var <allelic_variation_of_coding_mutation> -reg_mu <regulation_mutation_rate> -m_rate <migration_rate> -x_start <x_start> -y_start <y_start> -x_opt <x_optimum> -y_opt <y_optimum> -start_network <starting_network> -sel_var <variance_in_stabilizing_selection> -sel_covar <covariance_in_stabilizing_selection> -rec <recombination_rate> -theta <theta> -gamma <gamma> -model <model_A_or_B>

Example:
./simulation -npops 2 -i 1000 -gen 2 -mu 0.01 -mu_var 0.01 -reg_mu 0.01 -m_rate 0.1 -x_start 300 -y_start 300 -x_opt 300 -y_opt 300 -start_network 1111 -sel_var 100000 -sel_covar 0 -rec 0.35 -theta 300 -gamma 1 -model 'B'

Notes:
For random number generation C++11 is required


Developmental tools for troubleshooting memory leaks:
valgrind ./simulation -npops 2 -i 1000 -gen 2 -mu 0.01 -mu_var 0.01 -reg_mu 0.01 -m_rate 0.1 -x_start 300 -y_start 300 -x_opt 300 -y_opt 300 -start_network 1111 -sel_var 100000 -sel_covar 0 -rec 0.35 -theta 300 -gamma 1 -model 'B'


