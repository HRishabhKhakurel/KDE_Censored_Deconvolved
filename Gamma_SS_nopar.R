# Initial parameters
n = 2   # no. of repition
N = 10  # sample size
h_start = 0.1 # bandwidth sequence start
h_end = 0.2 # bandwidth sequence end

#sink(paste0('Gamma_SS_Output_',n,'_',N))

# Loading packages
print('-----------------------------------------------------------------')
print('Loading required R scripts')
source('KDE_decon_censored.R')
source('censoring_distribution.R')
source('eiv_data.R')
source('kernel_functions.R')
source('bandwidth_selection.R')
source('bandwidth_selection_parallel.R')
#source('multiResultClass.R')

print('-----------------------------------------------------------------')
print('Inputs')
print(paste0('Sample size: ',N))
print(paste0('No. of Repitions: ',n))
print(paste0('h start: ', h_start))
print(paste0('h_end: ', h_end))

# Setting up cluster
# print('-----------------------------------------------------------------')
# print('Setting Up Cluster')
# n.cores <- detectCores()
# print(paste0('no. of cores: ', n.cores))
# registerDoParallel(n.cores - 1)
# print('Cores Registered')

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
j = 1
while (j <= n){
        # generating data
        data = data.gen.gamma.SS.error(6,0.1,1,0.8,N,0.05)
        
        if(nrow(Cen_dist_KM(data$W_err,data$d2)) != N) {
        break 
        }
        else {

        # print(j)
        # finding error variance
        s <-  var(data$W_err) - var(data$W_err)*0.7
        
        # finding bandwidth
        h_opt_KM[j] <- lscv_cens(data,s,h_start, h_end,
                                      cens_dist = 'KM',kernel = 'normal')[,1]
        h_opt_BS[j] <- lscv_cens(data,s,h_start, h_end,
                                      cens_dist = 'BS',kernel = 'normal')[,1]
        
        
        # density estimation
        density_KM[,j] <- kde_cen_decon(data, h = h_opt_KM[j] ,s = s,
                                        cens_dist = 'KM',
                                        kernel = 'normal')[,2]
        density_BS[,j] <- kde_cen_decon(data, h = h_opt_BS[j],s = s,
                                        cens_dist = 'BS',
                                        kernel = 'normal')[,2]
        
        # calculating errors
        x <- seq(0,2,by = 0.025)
        f <- dgamma(x,shape = 6, scale = 0.1)
        final_mise_KM[j] <- trapz(x,(f- density_KM[,j])^2)
        final_mise_BS[j] <- trapz(x,(f- density_BS[,j])^2)
        j = j + 1
        }
}

print('-----------------------------------------------------------------')
print('Statistics:')
print(paste0('Mean ISE KM: ', mean(final_mise_KM, na.rm =  TRUE)))
print(paste0('Mean ISE BS: ', mean(final_mise_BS, na.rm = TRUE)))

# print('-----------------------------------------------------------------')
# print('Wrinting to CSV')
# write.csv(data.frame('KM' = h_opt_KM,'BS' = h_opt_BS),paste0('h_opt_gamma_SS_',n,'_',N,'.csv'))
# write.csv(data.frame('KM' = final_mise_KM,'BS' = final_mise_BS),paste0('final_mise_gamma_SS_',n,'_',N,'.csv'))
# write.csv(density_KM,paste0('density_KM_gamma_SS_',n,'_',N,'.csv'))
# write.csv(density_BS,paste0('density_BS_gamma_SS_',n,'_',N,'.csv'))

#sink()        