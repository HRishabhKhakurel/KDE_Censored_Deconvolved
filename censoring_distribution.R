# Censoring distribution estimator

# Kaplan Meier Based Estimator
library(survival)
Cen_dist_KM <- function(W,delta){
        ndelta <- 1 - delta
        s_fit <- survfit(Surv(W,ndelta)~1)
        G <- 1-s_fit$surv
        ind = which(G == 1)
        G[ind] = G[ind] - 0.0001
        df <- data.frame(sort(W),G)
        colnames(df) <- c('W_err','G')
        return(df)
}

# Blum and Sursula based estimator
Cen_dist_BS <- function(W,delta){
        ndelta <- 1 - delta
        N <- length(W)
        df <- data.frame(W,ndelta)
        S <- rep(1,N)
        Z <- df[order(df$W),]
        for (i in 1:N){
                if (W[i] <= Z[1,1]){
                        S[i] = 1
                }
                else if ((W[i] > Z$W[1]) & (W[i]<= Z$W[N])){
                        temp = Z[Z[,1] < W[i],]
                        for (j in 1:nrow(temp)){
                                S[i] = S[i]*((N-j+1)/(N-j+2))^(temp$ndelta[j])
                        }
                }
                else if (W[i] > Z[N,1]){
                        for (j in 1:nrow(temp)){
                                S[i] = S[i]*((N-j+1)/(N-j+2))^(df$ndelta[j])
                                }
                }
        }
        df1 <- data.frame(W,1-S)
        df2 <- df1[order(df1$W),]
        colnames(df2) <- c("W_err","S")
        return(df2)
}

