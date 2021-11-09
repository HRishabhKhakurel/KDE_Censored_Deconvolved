library(parallel)
library(doParallel)

source('KDE_decon_censored.R')
source('censoring_distribution.R')
source('eiv_data.R')
source('kernel_functions.R')
source('bandwidth_selection.R')
source('bandwidth_selection_parallel.R')

sink(file ='output.txt')
n.cores <- detectCores()
print(paste0('no. of cores: ', n.cores))
registerDoParallel(n.cores-1)
print('Cores Registered')

print('Generating data')
set.seed(6)
#generating data
data = data.gen.gamma.OS.error(6,0.1,1,0.8,50,0.05)
s <-  var(data$W_err) - var(data$W_err)*0.7
parallel_t <- system.time(h_opt <- lscv_cens_par(data,s,0.1,0.12,cens_dist = 'BS', kernel = 'double exponential'))
print('Parallel Run Time:')
print(parallel_t)
print('---------------------')
print('h_opt')
print(h_opt)
sink()