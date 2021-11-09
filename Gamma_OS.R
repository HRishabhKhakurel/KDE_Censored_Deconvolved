library(parallel)
library(doParallel)

# Initial parameters
n = 2   # no. of repition
N = 10  # sample size
h_start = 0.1 # bandwidth sequence start
h_end = 0.13 # bandwidth sequence end

# Loading packages
print('-----------------------------------------------------------------')
print('Loading required R scripts')
source('KDE_decon_censored.R')
source('censoring_distribution.R')
source('eiv_data.R')
source('kernel_functions.R')
source('bandwidth_selection.R')
source('bandwidth_selection_parallel.R')
source('multiResultClass.R')

print('-----------------------------------------------------------------')
print('Inputs')
sink(paste0('Gamma_OS_Output_',n,'_',N))
print(paste0('Sample size: ',N))
print(paste0('No. of Repitions: ',n))
print(paste0('h start: ', h_start))
print(paste0('h_end: ', h_end))

# Setting up cluster
print('-----------------------------------------------------------------')
print('Setting Up Cluster')
n.cores <- detectCores()
print(paste0('no. of cores: ', n.cores))
registerDoParallel(n.cores - 1)
print('Cores Registered')

# vectors to store results
print('-----------------------------------------------------------------')
print('Creating Required Vectors and Data Frames')
final_mise_KM = rep(NA,n)
final_mise_BS = rep(NA,n)
h_opt_KM = rep(NA,n)
h_opt_BS = rep(NA,n)
density_KM = data.frame(matrix(NA, nrow = 81, ncol = n))
density_BS = data.frame(matrix(NA, nrow = 81, ncol = n))

# Main Loop
print('-----------------------------------------------------------------')
print("Running Main Loop")
result <- foreach(j = 1:n) %dopar% {
        # generating data
        data = data.gen.gamma.OS.error(6,0.1,1,0.8,N,0.05)
        
        # skip_to_next <- FALSE
        # r <- tryCatch(kde_cen_decon(data, h = 0.09,s = s,
        #                             cens_dist = 'KM',
        #                             kernel = 'double exponential'),
        #               error = function(e) 
        #               { skip_to_next <<- TRUE})
        res <- multiResultClass()
        if(nrow(Cen_dist_KM(data$W_err,data$d2)) != N) {
                res$h_opt_KM = NA
                res$h_opt_BS = NA
                res$density_KM = rep(NA,N)
                res$density_BS = rep(NA,N)
                res$final_mise_KM = NA
                res$final_mise_BS = NA
        }
        else{
        # finding error variance
        s <-  var(data$W_err) - var(data$W_err)*0.7
        
        # finding bandwidth
        res$h_opt_KM <- lscv_cens_par(data,s,h_start, h_end,
                                     cens_dist = 'KM',kernel = 'double exponential')[,1]
        res$h_opt_BS <- lscv_cens_par(data,s,h_start, h_end,
                                     cens_dist = 'BS',kernel = 'double exponential')[,1]
        
        
        
        # density estimation
        res$density_KM <- kde_cen_decon(data, h = res$h_opt_KM ,s = s,
                                         cens_dist = 'KM',
                                         kernel = 'double exponential')[,2]
        res$density_BS <- kde_cen_decon(data, h = res$h_opt_BS,s = s,
                                            cens_dist = 'BS',
                                            kernel = 'double exponential')[,2]
        
        # res$density_KM <- result_KM$f_hat_scaled
        # res$density_BS <- result_BS$f_hat_scaled
        # calculating errors
        x <- seq(0,2,by = 0.025)
        f <- dgamma(x,shape = 6, scale = 0.1)
        res$final_mise_KM <- trapz(x,(f- res$density_KM)^2)
        res$final_mise_BS <- trapz(x,(f- res$density_BS)^2)
        }
#         return(c(h_opt_BS_temp, h_opt_KM_temp,
#                density_KM_temp,density_BS_temp,
#                final_mise_KM_temp, final_mise_KM_temp))
        return(res)
}
# mean_mise_KM = mean(final_mise_KM)
# mean_mise_BS = mean(final_mise_BS)

# Arranging Final Outputs
print('-----------------------------------------------------------------')
print('Arranging Final Outputs')
for (i in 1:n){
        h_opt_KM[i] = result[[i]]$h_opt_KM
        h_opt_BS[i] = result[[i]]$h_opt_BS
        final_mise_KM[i] = result[[i]]$final_mise_KM
        final_mise_BS[i] = result[[i]]$final_mise_BS
        density_KM[,i] = result[[i]]$density_KM
        density_BS[,i] = result[[i]]$density_BS
}

print('-----------------------------------------------------------------')
print('Statistics:')
print(paste0('Mean ISE KM: ', mean(final_mise_KM, na.rm =  TRUE)))
print(paste0('Mean ISE BS: ', mean(final_mise_KM, na.rm = TRUE)))

print('-----------------------------------------------------------------')
print('Wrinting to CSV')
write.csv(data.frame('KM' = h_opt_KM,'BS' = h_opt_BS),paste0('h_opt_gamma_OS_',n,'_',N,'.csv'))
write.csv(data.frame('KM' = final_mise_KM,'BS' = final_mise_BS),paste0('final_mise_gamma_OS_',n,'_',N,'.csv'))
write.csv(density_KM,paste0('density_KM_gamma_OS_',n,'_',N,'.csv'))
write.csv(density_BS,paste0('density_BS_gamma_OS_',n,'_',N,'.csv'))

sink()
