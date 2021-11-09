multiResultClass <- function(h_opt_KM= NULL,h_opt_BS= NULL,
                             density_KM = NULL,density_BS = NULL,
                             final_mise_KM = NULL,final_mise_BS = NULL)
{
        me <- list(
                h_opt_KM = h_opt_KM,
                h_opt_BS = h_opt_BS,
                density_KM = density_KM,
                density_BS = density_BS,
                final_mise_KM = final_mise_KM,
                final_mise_BS = final_mise_BS
        )
        
        ## Set the name for the class
        class(me) <- append(class(me),"multiResultClass")
        return(me)
}