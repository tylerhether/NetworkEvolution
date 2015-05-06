# Enumerate the parameter space
rm(list=ls())
library(plyr)
library(dplyr)


npops = 2
i = 1000
gen = 10000
mu = c(0.02)
mu_var = 0.01
reg_mu = c(0.002)
m_rate = c(0.0001)
x_start = 300
y_start = 300
x_opt = c(400)
y_opt = c(400)
start_network = '1111r'
sel_var = 10000 
sel_covar = 0 
rec = c(0.5)
theta = 300
gamma = 1
model = 'B' 
start_deviation = 10
selection_mode = 1
output_freq = 1000 
rep = 1:20

output_pheno_file = 'phenotypes' 
output_fitness_file = 'fitness'



dat <- expand.grid(
	npops=npops,
	i=i,
	gen=gen,
	mu=mu,
	mu_var=mu_var,
	reg_mu=reg_mu,
	m_rate=m_rate,
	x_start=x_start,
	y_start=y_start,
	x_opt=x_opt,
	y_opt=y_opt,
	start_network=start_network,
	sel_var=sel_var,
	sel_covar=sel_covar,
	rec=rec,
	theta=theta,
	gamma=gamma,
	model=model,
	start_deviation=start_deviation,
	selection_mode=selection_mode,
	output_freq=output_freq,
	rep=rep
)

dat$output_pheno_file <- paste(output_pheno_file, 1:dim(dat)[1], ".txt", sep="")
dat$output_fitness_file <- paste(output_fitness_file, 1:dim(dat)[1], ".txt", sep="")

# dat <- mutate(dat, output_pheno_file=paste(
# 		npops,	i,	gen,	mu,	mu_var,	reg_mu,	m_rate,	x_start,	
# 		y_start,	x_opt,	y_opt,	start_network,	sel_var,	sel_covar,	
# 		rec,	theta,	gamma,	model,	start_deviation,	selection_mode,	
# 		output_freq, output_pheno_file,sep="_"))
# dat <- mutate(dat, output_fitness_file=paste(
# 		npops,	i,	gen,	mu,	mu_var,	reg_mu,	m_rate,	x_start,	
# 		y_start,	x_opt,	y_opt,	start_network,	sel_var,	sel_covar,	
# 		rec,	theta,	gamma,	model,	start_deviation,	selection_mode,	
# 		output_freq, output_fitness_file,sep="_"))


write.table(file="./ParamSpace.txt", dat, row.names=FALSE, quote=FALSE, col.names=FALSE, sep=" ")


