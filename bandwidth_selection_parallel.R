library(parallel)
library(foreach)
library(doParallel)


# Least Squares Cross Validation for Right Censored Data (Marron and Padgett)
lscv_cens_par <- function(data,s,h_start,h_end,cens_dist,kernel){
  W <- data$W_err
  delta <- data$d2
  
  h <- seq(h_start, h_end, by = 0.01)
  s <-  var(data$W_err) - var(data$W_err)*0.7
  #error <- rep(NA, length(h))
  G <- switch (cens_dist,
               'BS' = Cen_dist_BS(W,delta),
               'KM' = Cen_dist_KM(W,delta)
  )
  temp_df <- merge(data,G,by = 'W_err')
  w <- (0.1)*exp(-0.1*W)
  temp_df2 <- data.frame(W,w)
  colnames(temp_df2) <- c('W_err','w')
  temp_df <- merge(temp_df,temp_df2, by ='W_err')
  error <- foreach (i=1:length(h), .combine =c, .packages = 'pracma') %dopar% {
    result1 <- kde_cen_decon_modified(data,h[i],s,cens_dist = cens_dist, kernel = kernel)
    f_hat_complete <- result1$f_hat_scaled
    f_hat_complete_squared_wt <- ((f_hat_complete)^2)*temp_df[,4]
    term1 <- trapz(result1$W,f_hat_complete_squared_wt)
    
    term2_vec <- foreach (j = 1:nrow(temp_df), .combine = c, .packages = 'pracma') %dopar% {
      (kde_cen_decon_LOO(data.frame('W_err' = temp_df$W_err,'d2' = temp_df$d2),
                             h[i],s,ind = j,
                             cens_dist = cens_dist,
                             kernel = kernel)[,2]*temp_df[j,2]*temp_df[j,4])/(1-temp_df[j,3])
    }
    term2 <- sum(term2_vec)
    term2 <- (2*term2)/(nrow(data))
    error <- term1 - term2
  }
  df <- data.frame(h,error)
  h_opt <- df[which(df$error == min(df$error)),]
  return(h_opt)
}