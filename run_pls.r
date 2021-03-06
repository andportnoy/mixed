library(feather)
library(optparse)
library(lme4pureR)
library(minqa)

option_list <- list(
    make_option("--data"),
    make_option("--formula"),
    make_option("--randomdata")
)

opt <- parse_args(OptionParser(option_list=option_list))

data <- read_feather(opt$data)
formula <- as.formula(readLines(opt$formula))

ll <- plsform(formula, data)

#f <- file(opt$randomdata, "rb")
#ll$theta <- readBin(f, "numeric", length(ll$theta))
#close(f)

devfun <- do.call(pls, ll)
lower <- ifelse(ll$theta,0,-Inf)
upper <- rep.int(Inf, length(ll$theta))

bobyqa(ll$theta, devfun, lower=lower, upper=upper)
