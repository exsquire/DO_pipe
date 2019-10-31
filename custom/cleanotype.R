options(warn = -1)
cat("Running Clean-o-Type...\n")
#Perform the sex intensity diagnostics
source("sexDiag.R")
#Process and clean the cross2 object
source("cleanCross.R")
#Generate scan inputs and intermediates for producing cleaned probability objects
source("computeScanInputs.R")
oldw <- getOption("warn")
options(warn = oldw)
