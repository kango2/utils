library(tidyverse)
##declare empty list
pval <- list()
##generate p-value of t-test between two random distributions with a mean of 0 and std-dev of 1, 10,000 times
for (i in 1:10000) {
    x <- rnorm(1000)
    y <- rnorm(1000)
    pval[i] <- t.test(x,y, alternative = "two.sided", var.equal = TRUE)$p.value
}

##convert to the vector
p.values <- unlist(pval)

##draw histogram of p-values
hist(p.values)

##check how many p-values are significant below 0.05
sum(p.values < .05)

##function for bonferroni correction

bonferroni.correction <- function(p, n){
    return(p * n)
}

##get bonferroni corrected p-values
corrected.p <- bonferroni.correction(p.values, length(p.values))

##check how many p-values are significant below 0.05 after correction
sum(unlist(corrected.p) < .05)

##draw corrected p-values of histogram

hist(corrected.p)
