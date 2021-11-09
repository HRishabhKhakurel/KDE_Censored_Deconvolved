# Deconvolved Kernel functions

# Double Exponential Error with Gaussian Kernel
double_exponential <- function(y,h,s){
        K <- (1/sqrt(2*pi))*exp(-y^2/2)*(1 - (s^2/(2*h^2))*(y^2 - 1))
        return(K)
}

# Normal Error from Erols Work
normal <- function(y,h,s){
        K <- (1/sqrt(2*pi*(1 - s^2/h^2)))*exp(-y^2/(2*(1 - s^2/h^2)))
        return(K)
}

# Triweight kernel with normal erro

triwght <- function(y,h,s){
        integrand <- function(t){ cos(rep(t,length(y))*y)*((1-rep(t,length(y))^2)^3)*exp(((s^2)*(rep(t,length(y))^2))/(2*h^2))}
        K <- (1/pi)*integrate(integrand,0,1)
}