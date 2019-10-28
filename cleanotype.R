options(warn = -1)
cat("Running Clean-o-Type...\n")

source("fstGen.R")
source("proc_cross.R")
source("compute.R")
oldw <- getOption("warn")
options(warn = oldw)
