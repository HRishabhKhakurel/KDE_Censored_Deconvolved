library(nimble)
data.gen.weibull.OS.error <- function(a1,b1,a2,b2,N, sigma){
        W <- rep(NA,N) # observed lifetimes
        W_err <- rep(NA,N)
        d1 <- rep(NA,N) # observed indicators
        d2 <- rep(NA,N)
        n = 1
        while(n <= N){
                # generating the censoring time
                c <- rweibull(1,shape = a2, scale = b2) # censoring time
                
                # Lifetime without errors
                x <- rweibull(1,shape = a1, scale = b1) # True lifetime
                W[n] <- min(x,c)
                d1[n] <- ifelse(c < x, 0, 1)
                
                # generating error from truncated double exponential
                e <- rdexp(1,location = 0, scale = sigma)
                if (e < 0){
                        next
                }
                z <- x + e # adding errors to true lifetime
                W_err[n] <- min(z,c)
                d2[n] <- ifelse(c < z, 0, 1)
                n = n + 1
        }
        error_data <- data.frame(W_err,d2)
        return(error_data)
}

data.gen.weibull.SS.error <- function(a1,b1,a2,b2,N,sigma){
        W <- rep(NA,N) # observed lifetimes
        W_err <- rep(NA,N) # observed lifetimes with error
        d1 <- rep(NA,N) # observed indicators
        d2 <- rep(NA,N) # observed indicators for W1_err
        n = 1
        while(n <= N){
                # generating the censoring time
                c <- rweibull(1,shape = a2, scale = b2) # censoring time
                
                # Lifetime without errors
                x <- rweibull(1,shape = a1, scale = b1) # True lifetime
                W[n] <- min(x,c)
                d1[n] <- ifelse(c < x, 0, 1)
                
                
                # generating error from truncated double exponential
                e <- rnorm(1,mean = 0, sd = sigma)
                if (e < 0){
                        next
                }
                z <- x + e # adding errors to true lifetime
                W_err[n] <- min(z,c)
                d2[n] <- ifelse(c < z, 0, 1)
                n = n + 1
        }
        error_data <- data.frame(W_err,d2)
        return(error_data)
}

data.gen.gamma.SS.error <- function(a1,b1,a2,b2,N, sigma){
        W <- rep(NA,N) # observed lifetimes
        W_err <- rep(NA,N) # observed lifetimes with error
        d1 <- rep(NA,N) # observed indicators
        d2 <- rep(NA,N) # observed indicators for W1_err
        n = 1
        while(n <= N){
                # generating the censoring time
                c <- rgamma(1,shape = a2, scale = b2) # censoring time
                
                # Lifetime without errors
                x <- rgamma(1,shape = a1, scale = b1) # True lifetime
                W[n] <- min(x,c)
                d1[n] <- ifelse(c < x, 0, 1)
                
                #sigma = (var(x) - var(x)*0.7)/0.7
                
                # generating error from truncated double exponential
                e <- rnorm(1,mean = 0, sd = sigma)
                if (e < 0){
                        next
                }
                z <- x + e # adding errors to true lifetime
                W_err[n] <- min(z,c)
                d2[n] <- ifelse(c < z, 0, 1)
                n = n + 1
        }
        error_data <- data.frame(W_err,d2)
        return(error_data)
}

data.gen.gamma.OS.error <- function(a1,b1,a2,b2,N,sigma){
        W <- rep(NA,N) # observed lifetimes
        W_err <- rep(NA,N) # observed lifetimes with error
        d1 <- rep(NA,N) # observed indicators
        d2 <- rep(NA,N) # observed indicators for W1_err
        n = 1
        while(n <= N){
                # generating the censoring time
                c <- rgamma(1,shape = a2, scale = b2) # censoring time
                
                # Lifetime without errors
                x <- rgamma(1,shape = a1, scale = b1) # True lifetime
                W[n] <- min(x,c)
                d1[n] <- ifelse(c < x, 0, 1)
                
                # generating error from truncated double exponential
                e <- rdexp(1,location = 0, scale = sigma)
                if (e < 0){
                        next
                }
                z <- x + e # adding errors to true lifetime
                W_err[n] <- min(z,c)
                d2[n] <- ifelse(c < z, 0, 1)
                n = n + 1
        }
        error_data <- data.frame(W_err,d2)
        return(error_data)
}