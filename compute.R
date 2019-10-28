#Intermediate object creation - to be blasted on the cluster
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


rm(list=ls())
use_ncores = 10
#---------------------------------------------
cat("\nCalling script: compute.R\n")
cat("\nDefault cores set to:", use_ncores)
#---------------------------------------------
pathIn <- list(
  CrossRaw = list.files(path = "../../inputs/",pattern = "cross", full.names = T, ignore.case = T),
  CrossPass1 = list.files(path = "../output/objects/",pattern = "pass1.rds", full.names = T),
  gprRaw = list.files(path = "../../outputs/", pattern = "gprRaw.rds", full.names = T),
  prefPath = list.files(path = "../../inputs/", pattern = "prefix", full.names = T, ignore.case = T)
)

#Check
if(!all(sapply(pathIn, function(x) length(x) > 0))){
  cat("Error: One or more input(s) could not be found.\n")
  print(sapply(pathIn, function(x) length(x) > 0))
}else{
  cat("\nFound the following inputs: \n")
  print(sapply(pathIn, function(x) length(x) > 0))
}

cross_Raw <- readRDS(pathIn[[1]])
cross_pass1 <- readRDS(pathIn[[2]])
#Generate GPR_Raw from crossRaw


cat("\nPulling gpr_Raw from DO_Pipe outputs folder.\nPlease wait...")
gpr_Raw <- readRDS(pathIn[[3]])
cat("Finished.\n")

prefix <- readRDS(pathIn[[4]])
#----------------------------------------------

#Generate e - most are largely negative - few are above 2, those 
#are what we want to scrutinize as they represent a large
#departure from what the multiparent gpr suggests
cat("\nCalculating genotyping error LOD scores...")

e <- calc_errorlod(cross_Raw, gpr_Raw, cores = use_ncores)
e <- do.call("cbind", e)
saveRDS(e, file = paste0("../output/objects/",prefix,"_e.rds"))

cat("Finished.\n")
cat("Source: Lincoln SE, Lander ES (1992)\nSystematic detection of errors in genetic linkage data.\nGenomics 14:604-610.")

#Examine samps with error LOD greater than 2
#Error_ind is the significant eLOD per sample
#describes # of markers with eLOD > 2 / total markers typed
#rowSums counts the number of TRUEs where e > 2
error_ind <- rowSums(e > 2)/n_typed(cross_Raw)*100

cat("\nGenerating Error LOD Sample Diagnostic...\n")
#Kick out the pre-clean errorlod
eLab <- paste0(names(error_ind)," (", round(error_ind,2),")")
eDF <- data.frame(samp = eLab,
                         Index = seq_along(eLab), 
                         val = round(error_ind,2))

e_thresh <- mean(eDF$val)+(3*sd(eDF$val))

#GGPLOT
eLOD <- ggplot(eDF, 
                   aes(x =Index, y = val))+
  geom_point(color = "firebrick3", 
             size = 3)+
  labs(title = "Error LOD - Sample Diagnostic",
       y = "Markers with Error LOD > 2 (%)",
       x = "Mouse")+ 
  geom_label_repel(data = eDF[eDF$val >= e_thresh,],
                   aes(label = samp),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   ylim = c(max(eDF$val), NA),
                   segment.color = 'grey50')+
  theme_bw()

eLOD %>% ggexport(filename=paste0("../output/plots/",prefix,"_eLODSamp_Diag.pdf"))
#--------------------------------------------
#Now pass1 clean 'e'
cat("\nFiltering 'g' and 'e' objects for Cross Pass1...")
g <- readRDS(paste0("../output/objects/",prefix,"_g.rds"))
e <- e[ind_ids(cross_pass1),]
g <- g[ind_ids(cross_pass1),]
cat("Finished.\n")

#------------------------------------
#Markers with errorLOD scores > 2
#For each marker, ask how many were > errorLOD 2 in all 
#subject. This gives the percent significant error for each marker
error_mar <- colSums(e > 2)/n_typed(cross_pass1, "marker")*100
saveRDS(error_mar, file = paste0("../output/objects/",prefix,"_error_mar.rds"))

cat("\nGenerating Error LOD Marker Diagnostic...\n")
errmarDF <- data.frame(samp = names(error_mar),
                   Index = seq_along(error_mar), 
                   val = error_mar) %>% arrange(desc(val))

eMarLOD <- ggplot(errmarDF, 
       aes(x =Index, y = val))+
  geom_point(color = "coral2", 
             size = 3)+
  labs(title = "Error LOD - Marker Diagnostic",
       y = "Error LOD > 2 Across Samples (%)",
       x = "Marker")+ 
  geom_label_repel(data = errmarDF[1:20,],
                   aes(label = samp),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  theme(axis.text.x = element_blank())+
  theme_bw()

eMarLOD %>% ggexport(filename = paste0("../output/plots/",prefix,"_eLODMark_Diag.pdf"))
#------------------------------------
cat("\nGenerating Marker Tri-plot Diagnostic...")

fgn <- readRDS(paste0("../output/objects/",prefix,"_fgn.rds"))
gf_mar <- t(apply(g, 2, function(a) table(factor(a, 1:3))/sum(a != 0)))
gn_mar <- t(apply(g, 2, function(a) table(factor(a, 1:3))))


pdf(paste0("../output/plots/",prefix,"_MarkTriPlot.pdf"))
par(mfrow=c(2,2), mar=c(0.6, 0.6, 2.6, 0.6))
for(i in 1:4) {
  triplot(c("AA", "AB", "BB"), main=paste0("MAF = ", i, "/8"))
  z <- gf_mar[fgn==i,]
  z <- z[rowSums(is.na(z)) < 3,]
  tripoints(z, pch=21, bg="gray80", cex=0.6)
  tripoints(c((1-i/8)^2, 2*i/8*(1-i/8), (i/8)^2), pch=21, bg="violetred")
}
dev.off()
cat("Finished.\n")
#------------------------------------

cat("\nGenerating and Exporting Maximal Marginal Probability Genotypes...\n")
m <- maxmarg(gpr_Raw, minprob = 0.5, cores = use_ncores)
saveRDS(m, file = paste0("../output/objects/",prefix,"_m.rds"))
cat("Finished.\n")


cat("\nGenerating Crossover Diagnostic...\n")
nxo <- count_xo(m, cores = 4)
totxo <- rowSums(nxo)


xolab <- paste0(names(totxo), " (", totxo,")")
xoDF <- data.frame(samp = xolab,
                         Index = seq(length(xolab)), 
                         val = totxo)

xo_thresh <- mean(xoDF$val)+(3*sd(xoDF$val))

#GGPLOT
xoPlot <- ggplot(xoDF, 
                   aes(x =Index, y = val))+
  geom_point(color = "springgreen3", 
             size = 3)+
  labs(title = "Crossover Diagnostic",
       y = "Crossovers (n)",
       x = "Mouse")+ 
  geom_label_repel(data = xoDF[xoDF$val >= xo_thresh,],
                   aes(label = samp),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   ylim = c(max(xoDF$val), NA),
                   segment.color = 'grey50')+
  theme_bw()

xoPlot %>% ggexport(filename = paste0("../output/plots/",prefix,"_xoPlot.pdf")) 

#-----------------------------------
#generate cross_Clean 
cat("\nFiltering 'Cross_Pass1' for markers with > 5% Error Rate\n")
cross_Clean <- drop_markers(cross_pass1, names(error_mar)[error_mar > 5])

cat("Exporting 'Cross_Clean' object...")
saveRDS(cross_Clean,paste0("../../outputs/",prefix,"_cross_Clean.rds"))
cat("Finished.\n")

#Make gpr_CLEAN

cat("\nCalculating cleaned genotype probabilities.\nPlease wait...")
gpr_Clean <- calc_genoprob(cross_Clean, cores = use_ncores)
saveRDS(gpr_Clean, file = paste0("../output/objects/",prefix,"_gpr_Clean.rds"))
cat("Finished.\n")

cat("\nConverting to Allele Probabilities.\nPlease wait...")
apr_Clean <- genoprob_to_alleleprob(gpr_Clean, cores = use_ncores)
saveRDS(apr_Clean, file = paste0("../../outputs/",prefix,"_apr_Clean.rds"))
cat("Finished.\n")

cat("\nCalculating kinship.\nPlease wait...")
kLOCO_Clean <- calc_kinship(apr_Clean, type = "loco", cores = use_ncores)
saveRDS(kLOCO_Clean, file = paste0("../../outputs/",prefix,"_kLOCO_Clean.rds"))
cat("Finished.\n")

rm(apr_Clean)

cat("\nGenerating and Exporting Multiparent-predicted Genotypes...")
snpg <- predict_snpgeno(cross_Clean, m, cores = use_ncores)
snpg <- do.call("cbind", snpg)
saveRDS(snpg, file = paste0("../output/objects/",prefix,"_snpg.rds"))
cat("Finished.\n")


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


invisible(gc())
cat("\ncompute.R: Status: Complete")





























