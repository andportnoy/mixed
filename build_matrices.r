suppressPackageStartupMessages(library(lme4))
library(feather)
library(optparse)

option_list <- list(
    make_option("--data"),
    make_option("--formula"),
    make_option("--X"),
    make_option("--Z"),
    make_option("--Lambdatx"),
    make_option("--Lambdati"),
    make_option("--Lambdatp")
)
opt <- parse_args(OptionParser(option_list=option_list))

data <- read_feather(opt$data)
formula <- readLines(opt$formula)

lf <- lFormula(formula, data)

file <- file(opt$X, "wb")
writeBin(c(t(lf$X)), file)
close(file)

file <- file(opt$Z, "wb")
writeBin(c(t(as.matrix(lf$reTrms$Zt))), file)
close(file)

Lambdat <- lf$reTrms$Lambdat
file <- file(opt$Lambdatx, "wb")
writeBin(Lambdat@x, file)
close(file)
file <- file(opt$Lambdati, "wb")
writeBin(Lambdat@i, file)
close(file)
file <- file(opt$Lambdatp, "wb")
writeBin(Lambdat@p, file)
close(file)
