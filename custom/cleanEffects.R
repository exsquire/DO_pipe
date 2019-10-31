library(qtl2)
library(dplyr)
library(tibble)
library(broman)
library(ggpubr)
library(ggplot2)
library(ggrepel)
library(qtlcharts)
library(gridExtra)
library(data.table)

#Separate from pipeline and return as needed. This file should aggregate
#the different methods of comparing the raw and cleaned data. 

#Begin work by reading in the inputs
#Things to include - post XO and post error LOD


cat("\nAligning gpr_Raw and gpr_Clean...\n")
#Make gprs flush
gpr_Clean <- gpr_Clean[ind_ids(cross_Clean),]
gpr_Raw <- gpr_Raw[ind_ids(cross_Clean),]

for(i in seq_along(gpr_Raw)){
  message(i, " of ", length(gpr_Raw))
  tmpRaw <- colnames(gpr_Raw[[i]][1,,])
  tmpClean <- colnames(gpr_Clean[[i]][1,,])
  keep <- intersect(tmpRaw, tmpClean)
  gpr_Clean[[i]] <- gpr_Clean[[i]][,,keep]
  gpr_Raw[[i]] <- gpr_Raw[[i]][,,keep]
}
cat("Finished.\n")
cat("\nCalculating Genotype Probability Differences...\n")
#Find areas of largest disagreement
prdiff <- vector("list",length(gpr_Raw))
for(i in seq_along(prdiff)){
  message(i, " of ", length(prdiff))
  prdiff[[i]] <- apply(abs(gpr_Raw[[i]] - gpr_Clean[[i]]),c(1,3),sum)
}
names(prdiff) <- names(gpr_Raw)
cat("Finished.\n")



#----------------------------------------
#Differences take values between 0 and 2. 2 == completely different
#Show # of markers x samples per chromosome with differences > 1.5 (arbitrary, ~75% different)
cat("\nExporting Genotype Probability Difference Summaries...")
checkDiff_mark <- sapply(prdiff, function(d) sum(d > 1.5))
#Show # of samples per chromosome with at least 5 markers with an absolute difference > 1.5
checkDiff_samp <- sapply(prdiff, function(d) sum(rowSums(d > 1.5) > 5))
saveRDS(checkDiff_mark, file = paste0("../output/objects/",prefix,"_checkDiff_mark.rds"))
saveRDS(checkDiff_samp, file = paste0("../output/objects/",prefix,"_checkDiff_samp.rds"))
cat("Finished.\n")
#Perform Marker and Sample Specific Diagnostics in another script
#Pre-clean genotypes are in blue, post-clean is in red
#Pre and Post color scales are combined so agreement is in purple
#How to read: Look for blue, that's what HMM thought it was before
#cleaning. Scan vertically from blue to red, that is how HMM
#"Changed its mind" based on the cleaned data. 

#Roll up a threshold of pr diffs per chr
#The samples included have summed pr diffs greater than
#3sd from the mean of all samples in that chromosome
#Make a df showing samples and chrs past this threshold

compDF <- data.frame(samp = character(),
                     chr = integer(),
                     err = numeric(), 
                     stringsAsFactors = F,
                     row.names = NULL)

#Percent Error represent the number of markers with more than 
#75% disagreement
for(i in seq_along(prdiff)){
  tmp     <- rowSums(prdiff[[i]])
  thresh  <- mean(tmp) + 3*sd(tmp)
  check   <- names(which(tmp > thresh))  
  percErr <- numeric()
  for(j in seq_along(check)){
    err <- table(prdiff[[i]][check[[j]],] > 1.5)
    if(is.na(err[2])){
      percErr <- c(percErr, 0)
    }else{
      percErr <- c(percErr, round(err[2]/(sum(err)),4)*100)
    }
    
  }
  compDF  <- rbind(compDF, 
                   data.frame(samp = check,
                              chr = rep(i,length(check)),
                              err = percErr,
                              stringsAsFactors = F,
                              row.names = NULL))
  
}

#Subset compDF to >= 0.1 percent disagreement
compDF <- compDF[compDF$err > 0.1,]

cat("\nCreating plot output sub directory...")

cat("\nPlotting",nrow(compDF),"genotype comparisons...\n")
dir.create(path = "../output/plots/CleanEffDiag", showWarnings = F)
pmap <- cross_Clean$pmap
for(i in 1:nrow(compDF)){
  message(i, " of ", nrow(compDF))
  id <- compDF[i,"samp"]
  chr <- as.character(compDF[i,"chr"])
  if(chr == "20"){
    chr <- "X"
  }
  err <- compDF[i,"err"]
  desig <- paste0(prefix,"_Chr_",chr,
                  "_",id)
  pdf(file = paste0("../output/plots/CleanEffDiag/",desig,".pdf"))
  if(isdebugged(plot_genoprobcomp)){
    undebug(plot_genoprobcomp)
  }
  
  qtl2::plot_genoprobcomp(probs1 = gpr_Raw,
                          probs2 = gpr_Clean,
                          map = pmap,
                          ind = id, chr = chr,
                          threshold = 0.25,
                          main = paste0(desig," - Error: ",err,"%")) 
  dev.off()
}
cat("Finished.\n")
