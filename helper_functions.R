################################################################################
# Helper functions
################################################################################


getX <- function(ntpts, tint){
    # Intercept and intervention term
    cbind(1, cumsum(1:ntpts > tint) / 100)
}
