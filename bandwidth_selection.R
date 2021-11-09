# Bandwidth selection 

# Least Squares Cross Validation

lscv <- function(data,h_start,h_end,cens_dist,kernel){
        
        # sequence of bandwidth
        h <- seq(h_start, h_end, by = 0.01)
        
        # error variance 
        s <-  var(data$W_err) - var(data$W_err)*0.7
        
        # vector to store error
        error <- rep(NA, length(h))
        
        for (i in 1:length(h)){
                print(paste0('bandwidth no. ', i))
                # Term 1
                # Full KDE result
                result1 <- kde_cen_decon_modified(data,h[i],s,cens_dist = cens_dist, kernel = kernel)
                f_hat_complete <- result1$f_hat_scaled
                f_hat_complete_squared <- (f_hat_complete)^2 # probability density estimate squared
                term1 <- trapz(result1$W,f_hat_complete_squared)
                
                # Term 2
                term2 <- 0
                for (j in 1:nrow(data)){
                        print(paste0('LOOCV: ', j))
                        r <- kde_cen_decon_LOO(data,h[i],s,ind = j,cens_dist = cens_dist, kernel = kernel)
                        term2 <- term2 + r[,2]
                }
                term2 <- (2*term2)/(nrow(data))
                error[i] <- term1 - term2
        }
        return(data.frame(h,error))
}


# Least Squares Cross Validation for Right Censored Data (Marron and Padgett)
lscv_cens <- function(data,s,h_start,h_end,cens_dist,kernel){
        W <- data$W_err
        delta <- data$d2
        
        h <- seq(h_start, h_end, by = 0.01)
        s <-  var(data$W_err) - var(data$W_err)*0.7
        error <- rep(NA, length(h))
        G <- switch (cens_dist,
                'BS' = Cen_dist_BS(W,delta),
                'KM' = Cen_dist_KM(W,delta)
        )
        temp_df <- merge(data,G,by = 'W_err')
        w <- (0.1)*exp(-0.1*W)
        temp_df2 <- data.frame(W,w)
        colnames(temp_df2) <- c('W_err','w')
        temp_df <- merge(temp_df,temp_df2, by ='W_err')
        for (i in 1:length(h)){
                # print(paste0('bandwidth no. ', i))
                result1 <- kde_cen_decon_modified(data,h[i],s,cens_dist = cens_dist, kernel = kernel)
                f_hat_complete <- result1$f_hat_scaled
                f_hat_complete_squared_wt <- ((f_hat_complete)^2)*temp_df[,4]
                term1 <- trapz(result1$W,f_hat_complete_squared_wt)
                
                term2 <- 0
                for (j in 1:nrow(temp_df)){
                        # print(paste0('LOOCV: ', j))
                        r <- kde_cen_decon_LOO(data.frame('W_err' = temp_df$W_err,'d2' = temp_df$d2),
                                               h[i],s,ind = j,
                                               cens_dist = cens_dist,
                                               kernel = kernel)
                        term2 <- term2 + (r[,2]*temp_df[j,2]*temp_df[j,4])/(1-temp_df[j,3])
                        }
                term2 <- (2*term2)/(nrow(data))
                error[i] <- term1 - term2
                # print(error)
        }
        df <- data.frame(h,error)
        h_opt <- df[which(df$error == min(df$error)),]
        return(h_opt)
}

# likelihood crossvalidation
llcv <- function(data,h_start,h_end,cens_dist,kernel){
        h <- seq(h_start, h_end, by = 0.01)
        s <-  var(data$W_err) - var(data$W_err)*0.7
        error <- rep(NA, length(h))
        
        for (i in 1:length(h)){
                print(paste0('bandwidth no. ', i))
                sum <- 0
                for (j in 1:nrow(data)){
                        print(paste0('LOOCV: ', j))
                        r <- kde_cen_decon_LOO(data,h[i],s,ind = j,cens_dist = cens_dist, kernel = kernel)
                        sum <- sum + log(r[,2])
                }
                error[i] <- -sum/nrow(data)
        }
        df <- data.frame(h,error)
        h_opt <- df[which(df$error == min(df$error)),]
        return(h_opt)
}


# term2 <- 0
# for (j in 1:nrow(data)){
#         r <- kde_cen_decon_LOO(data,0.11,s,ind = j,cens_dist = 'KM', kernel = 'double exponential')                 
#         print(paste0('r: ',r[,2]))
#         term2 <- term2 + (r[,2]*temp_df[j,2]*temp_df[j,4])/(1-temp_df[j,3])
# }
# 
# term2 <- 0
# for (j in 1:nrow(temp_df)){
#         print(paste0('LOOCV: ', j))
#         r <- kde_cen_decon_LOO(data.frame('W_err' = temp_df$W_err,'d2' = temp_df$d2),
#                                0.11,s,ind = j,
#                                cens_dist = 'KM',
#                                kernel = 'double exponential')
#         print(paste0('r: ',r[,2]))
#         term2 <- term2 + (r[,2]*temp_df[j,2]*temp_df[j,4])/(1-temp_df[j,3])
#         print(paste0('term2: ',term2))
# }