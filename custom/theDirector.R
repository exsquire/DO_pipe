#permAlloc: Tool for predicting time and memory usage of sbatch jobs
#Load requisite scan inputs
library(qtl2)

pathIn <- list(
  aprPath = list.files(path = "../outputs/",pattern = "apr_Clean", full.names = T, ignore.case = T),
  crossPath = list.files(path = "../outputs/",pattern = "cross_Clean", full.names = T, ignore.case = T),
  kinPath = list.files(path = "../outputs/", pattern = "kLOCO_Clean", full.names = T),
  covPath = list.files(path = "../inputs/", pattern = "covar", full.names = T, ignore.case = T),
  prefPath = list.files(path = "../inputs/", pattern = "prefix", full.names = T, ignore.case = T)
)

#Check
if(!all(sapply(pathIn, function(x) length(x) > 0))){
  cat("Error: One or more input(s) could not be found.\n")
  print(sapply(pathIn, function(x) length(x) > 0))
}else{
  cat("\nFound the following inputs: \n")
  print(sapply(pathIn, function(x) length(x) > 0))
}

cat("Loading Inputs, please wait...\n")
apr <- readRDS(pathIn[[1]])
cross <- readRDS(pathIn[[2]])
kLOCO <- readRDS(pathIn[[3]])
covar <- readRDS(pathIn[[4]])
prefix <- readRDS(pathIn[[5]])

ncores = 4

#The memory allocation should be mem_used() + max memory used by function
nPerm <- c(1, 5, 10, 20)

timeList <- list()
memList <- list()
cat("Estimating memory and time allocation...\n")
cat("Takes around 15-20 min, go do something else.\n")
for(i in 1:4){
  cat(nPerm[i], "permutation(s)\n")
  invisible(gc(reset = T))
  t <- system.time(
    scan1perm(apr, cross$pheno[,2,drop = F],
              kinship = kLOCO, 
              addcovar = covar,
              n_perm = nPerm[i],
              cores = ncores)
  )
  memList[[i]]<- sum(gc()[,6])
  timeList[[i]]<- as.numeric(t["elapsed"])
}
cat("Done.\n")

#Fit a predictive model that regresses time on number of permutations
res <- data.frame(perms = nPerm,
                  maxMem = sapply(memList, function(x) return(x)),
                  time = sapply(timeList, function(x) return(x)))
fit <- lm(time ~ perms ,data = res)

#Determine the number of permutations that finish around 10 minutes
#and round up to nearest value evenly divisible by 1000
lim <- 0
tryPerm <- 1
while(lim < 10){
  lim <- predict(fit, data.frame(perms = tryPerm))/60
  tryPerm <- tryPerm + 1
}

while(1000 %% tryPerm != 0){
  tryPerm <- tryPerm + 1
}

#Estimate the time for that number of permutations = 1 job
#Time is not a limiting factor and can be unstable on the FARM, so overestimate by 100% to account for I/O
timEst <- predict(fit, data.frame(perms = tryPerm)) * 2
#Overestimate the maxMem allocation by 10% and round up to nearest gig. 
memEst <- ceiling((max(res$maxMem) * 1.1)/1000)*1000

#Use memEst*4cores/maxCoreMem to determine the number of cores to request
needCores <- ceiling(memEst*4/2500)

#Estimate the number of parallelisms based on estimate memory use and max group limits
parJobs <- floor(250000/(needCores * 2500))
totJobs <- ncol(cross$pheno)*(1000/tryPerm)
jobMin <- ceiling(timEst/60)
runTime <- ceiling(totJobs*jobMin/parJobs)


cat("Assuming [" ,totJobs, "] jobs, run [",
    parJobs, "] at a time, \nwith each job lasting [", 
    jobMin, "] minutes, \nthe total array will complete in ~ [", 
    paste0(floor(runTime/60), " hr ",
           floor((runTime/60 - floor(runTime/60))*60), " min"), "] minutes.\n")

timeOut <- paste0(formatC(ceiling(timEst/60/60),width = 2, flag = "0"),":00:00")

permPipe <- function(x, numCols = 1, numPerms){
  #Submit phenotype matrix - ask for the col and perm parameters
  #define cuts 
  starts <- rep(seq(1,ncol(x), by = numCols), each = 1000/tryPerm)
  stops <-  c(starts[1:(length(starts)-(1000/tryPerm))] + (numCols-1),
              rep(ncol(x),1000/tryPerm))
  #No need for seed if the permutations aren't split
  aID <- seq(1, length(starts))
  
  #build control object
  ctrl <- data.frame(aID = aID,
                     start = starts,
                     stop = stops,
                     n_perm = numPerms)
  #Output messages
  cat("\nControl file will generate [", nrow(ctrl), "] jobs.")
  cat("\nA single job consists of [", numCols, "] phenotype(s),\nand [",numPerms,"] permutations.\n")
  return(ctrl)
}

#Genrate controller file for sbatch
ctrl <- permPipe(cross$pheno, 1, tryPerm)
saveRDS(ctrl, file = paste0(prefix, "_controller.rds"))

#Create sbatch argument file
sink("batchArgs.txt")
cat(tryPerm, 
    ncol(cross$pheno) * (1000/tryPerm),
    timeOut,
    needCores)
sink()

