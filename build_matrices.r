suppressPackageStartupMessages(library(lme4))
library(feather)
library(optparse)

option_list <- list(
    make_option("--data"),
    make_option("--formula"),
    make_option("--X"),
    make_option("--Z"),
    make_option("--L")
)
opt <- parse_args(OptionParser(option_list=option_list))

data <- read_feather(opt$data)
formula <- readLines(opt$formula)

lf = lFormula(formula, data)

file <- file(opt$X, "wb")
writeBin(c(t(lf$X)), file)
close(file)

file <- file(opt$Z, "wb")
writeBin(c(t(as.matrix(lf$reTrms$Zt))), file)
close(file)

file <- file(opt$L, "wb")
writeBin(c(as.matrix(lf$reTrms$Lambdat)), file)
close(file)
