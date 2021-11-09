# Deconvolved Density Estimator for Right Censored Data
library(pracma)

kde_cen_decon <- function(data,h,s,
                          kernel = c('double exponential','normal','triwght'),
                          cens_dist = c('KM','BS'))
        {
        
        # Extracting some required data
        W <- data[,1]
        delta <- data[,2]
        
        # length of data
        N <- length(W)
        
        # vectors to store required resutls
        K <- vector() # vector to store values from kernel functions
        M <- vector() # vector to store values from inverse censoring weights
        f_hat <- vector() # vector to store estimated densities
        
        # estimating censoring distribution
        # KM uses Kaplan Meier
        # BS uses Blum and Sursula
        G <- switch(cens_dist,
                    'KM' = Cen_dist_KM(W,delta),
                    'BS' = Cen_dist_BS(W,delta))
        temp_df <- merge(data,G,by = 'W_err')
        for (i in 1:N){
                M[i] <- temp_df[i,2]/(1-temp_df[i,3])
        }
        M_mat <- matrix(M,nrow = length(M),ncol = 1)
        
        
        #creating some other required vectors
        x <- seq(0,2, by =0.025)
        
        for (i in 1:length(x)){
                y <- vector()
                for (j in 1:N){
                        y[j] <- (x[i] - temp_df[j,1])/h       #(x-X[i])/h
                }
                K <- switch(kernel,
                            'double exponential' = double_exponential(y,h,s),
                            'normal' = normal(y,h,s),
                            'triwght' = triwght(y,h,s))
                K[K < 0] <- 0
                K_mat <- matrix(K,nrow = 1, ncol = length(K))
                f_hat[i] <- ((1/(N*h))*(K_mat %*% M_mat))[1,1]
        }
        
        # Standardizing the integral over the density to be equal to 1
        i <- sum(f_hat * 0.05)
        f_hat_scaled <- f_hat/i
        
        return(data.frame(x,f_hat_scaled))
}

kde_cen_decon_modified <- function(data,h,s,
                          kernel = c('double exponential','normal','triwght'),
                          cens_dist = c('KM','BS'))
{
        
        # Extracting some required data
        W <- data[,1]
        delta <- data[,2]
        
        # length of data
        N <- length(W)
        
        # vectors to store required resutls
        K <- vector() # vector to store values from kernel functions
        M <- vector() # vector to store values from inverse censoring weights
        f_hat <- vector() # vector to store estimated densities
        
        # estimating censoring distribution
        # KM uses Kaplan Meier
        # BS uses Blum and Sursula
        G <- switch(cens_dist,
                    'KM' = Cen_dist_KM(W,delta),
                    'BS' = Cen_dist_BS(W,delta))
        temp_df <- merge(data,G,by = 'W_err')
        for (i in 1:N){
                M[i] <- temp_df[i,2]/(1-temp_df[i,3])
        }
        M_mat <- matrix(M,nrow = length(M),ncol = 1)
        
        #creating some other required vectors
        #x <- seq(floor(min(W)),ceiling(max(W)), by =0.05)
        #W_sort <- sort(W)
        for (i in 1:length(temp_df[,1])){
                y <- vector()
                for (j in 1:N){
                        y[j] <- (temp_df[i,1] - temp_df[j,1])/h       #(x-X[i])/h
                }
                K <- switch(kernel,
                            'double exponential' = double_exponential(y,h,s),
                            'normal' = normal(y,h,s),
                            'triwght' = triwght(y,h,s))
                K[K < 0] <- 0
                K_mat <- matrix(K,nrow = 1, ncol = length(K))
                f_hat[i] <- ((1/(N*h))*(K_mat %*% M_mat))[1,1]
        }
        
        # Standardizing the integral over the density to be equal to 1
        i <- trapz(temp_df[,1],f_hat)
        f_hat_scaled <- f_hat/i
        final_df <- data.frame(temp_df[,1],f_hat_scaled)
        colnames(final_df) = c('W','f_hat_scaled')
        return(final_df)
}

kde_cen_decon_LOO <- function(data,h,s,ind,
                            kernel = c('double exponential','normal','triwght'),
                            cens_dist = c('KM','BS'))
{
        
        # Extracting some required data
        W <- data[,1]
        delta <- data[,2]
        
        # leave out the data indexed "ind"
        # W_loo <- W[-ind]
        # delta_loo <- W[-ind]
        
        # length of data
        # N <- length(W_loo)
        
        # vectors to store required resutls
        K <- vector() # vector to store values from kernel functions
        M <- vector() # vector to store values from inverse censoring weights
        f_hat <- vector() # vector to store estimated densities
        
        data_loo <- data[-ind,]
        N <- nrow(data_loo)
        # estimating censoring distribution
        # KM uses Kaplan Meier
        # BS uses Blum and Sursula
        # print(paste0('time: ',length(data_loo[,1])))
        # print(paste0('status: ',length(data_loo[,2])))
        #print(Cen_dist_KM(data_loo[,1],data_loo[,2]))
        G <- switch(cens_dist,
                    'KM' = Cen_dist_KM(data_loo[,1],data_loo[,2]),
                    'BS' = Cen_dist_BS(data_loo[,1],data_loo[,2]))
        temp_df <- merge(data_loo,G,by = 'W_err')
        for (i in 1:N){
                M[i] <- temp_df[i,2]/(1-temp_df[i,3])
        }
        M_mat <- matrix(M,nrow = length(M),ncol = 1)
        
        #creating some other required vectors
        #x <- seq(floor(min(W)),ceiling(max(W)), by =0.05)
        W_sort <- sort(W)
        for (i in 1:length(W_sort)){
                y <- vector()
                for (j in 1:N){
                        y[j] <- (W_sort[i] - temp_df[j,1])/h       #(x-X[i])/h
                }
                K <- switch(kernel,
                            'double exponential' = double_exponential(y,h,s),
                            'normal' = normal(y,h,s),
                            'triwght' = triwght(y,h,s))
                K[K < 0] <- 0
                K_mat <- matrix(K,nrow = 1, ncol = length(K))
                f_hat[i] <- ((1/(N*h))*(K_mat %*% M_mat))[1,1]
        }
        
        # Standardizing the integral over the density to be equal to 1
        i <- trapz(W_sort,f_hat)
        f_hat_scaled <- f_hat/i
        df <- data.frame(W_sort,f_hat_scaled)
        ind2 <- which(round(df$W_sort,6) == round(W[ind],6))
        return(df[ind2,])
}