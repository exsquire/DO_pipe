#Functions
missingNo <- function(pathIn, numJobs){
  actual <- as.numeric(gsub("[^0-9]*","",
                            gsub("^.*_perm_","",list.files(pathIn))))
  expected <- 1:numJobs
  return(expected[!expected %in% actual])
}

quiltR <- function(pathIn){
  library(utils)
  files <- list.files(pathIn, full.names = T)
  #Init output list
  outList <- list()
  #Method: Pull in any tmp file, ask it for it's column names, loop through the column names
  #ask if the column name exists as a slot in the list. If it does, append, if not, add.
  pb <- txtProgressBar(min = 0, max = length(files), style = 3)
  for(i in seq_along(files)){
    #Pull the perm output
    tmp <- readRDS(files[i])
    #seq along the columns
    for(j in seq_along(colnames(tmp))){
      if(!colnames(tmp)[j] %in% names(outList)){
        outList[[colnames(tmp)[j]]] <- tmp[[j]]
      }else{
        #concat output to slot data and set
        outList[[colnames(tmp)[j]]] <- c(outList[[colnames(tmp)[j]]],tmp[[j]])
      }
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  out <- do.call("cbind", outList)
}

#quiltR for DO_pipe
pathIn <- "../outputs"
prefix <- readRDS("../../inputs/prefix.rds")

#Load in batchArgs.txt
numJobs <- as.numeric(unlist(strsplit(readLines("../inputs/batchArgs.txt", warn = F), split = " "))[2])



failed <- missingNo(pathIn, numJobs)

if(length(failed) != 0){
  stop(paste0("Uh oh! The following jobs failed:\n", failed,"\nRe-run them and try again.\n"))
}else{
  permMat <- quiltR(pathIn)
  saveRDS(permMat, file = paste0("../../outputs/",prefix,"_fullPerm.rds"))
  permThresh <- apply(permMat, 2, quantile, probs = c(0.63, 0.95, 0.99))
  saveRDS(permThresh, file = paste0("../../outputs/",prefix, "_permThresh.rds"))          
}

