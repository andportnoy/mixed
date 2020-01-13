library(feather)
library(optparse)
library(lme4pureR)

option_list <- list(
    make_option("--data"),
    make_option("--formula")
)

opt <- parse_args(OptionParser(option_list=option_list))

data <- read_feather(opt$data)
formula <- as.formula(readLines(opt$formula))

ll <- plsform(formula, data)
#print(ll)

thetafun <- do.call(pls, ll)
#print(thetafun(ll$theta))

#thetafun2 <- pls(ll$X, ll$y, ll$Zt, ll$Lambdat, ll$thfun)
#print(thetafun2(ll$theta))

