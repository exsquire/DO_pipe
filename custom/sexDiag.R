library(fst)
library(qtl2)
library(broman)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(qtlcharts)
library(data.table)

rm(list=ls())
cat("\nCalling script: sexDiag.R\n")

#Run script from folder with the 3 input files 
#Form a list of filepaths where [[1]] == cross2, 
#[[2]] == codes, and [[3]] == Final Report
pathIn <- list(
  crossPath = list.files("../../inputs/", pattern = "cross", full.names = T, ignore.case = T),
  codePath  = list.files("../../inputs/", pattern = "allelecodes", full.names = T, ignore.case = T),
  reportPath= list.files("../../inputs", pattern = "finrep", full.names = T, ignore.case = T),
  prefPath = list.files("../../inputs/", pattern = "prefix", full.names = T, ignore.case = T)
)


if(!all(sapply(pathIn, function(x) length(x) > 0))){
  cat("Error: One or more inputs could not be found.\n")
  print(sapply(pathIn, function(x) length(x) > 0))
}else{
  cat("\nFound the following inputs: \n")
  print(sapply(pathIn, function(x) length(x) > 0))
}

#Load cross
cross <- readRDS(pathIn[[1]])

#Get sex status - both or single
if(length(unique(cross$covar$sex)) == 1){
  sexStat <- "single"
}else{
  sexStat <- "both"
}

#Load codes
codes <- read.csv(pathIn[[2]], skip = 3, stringsAsFactors = F) 
#Load a final report
cat("\nReading in Final Report. Please wait...\n")
finrep <- fread(pathIn[[3]], data.table = FALSE)
prefix <- readRDS(pathIn[[4]])

cat("Pulling Intensities...\n")

# create matrices that are snps x samples
snps <- unique(finrep[,"SNP Name"])
samples <- unique(finrep[,"Sample ID"])
X <- Y <- matrix(ncol=length(samples), nrow=length(snps))
dimnames(X) <- dimnames(Y) <- list(snps, samples)
for(i in seq(along=samples)) {
  message(i, " of ", length(samples))
  tmp <- finrep[finrep[,"Sample ID"]==samples[i],]
  X[,samples[i]] <- tmp[,"X"]
  Y[,samples[i]] <- tmp[,"Y"]
}
cat("\nFinished. \n")
# bring together in one matrix
result <- cbind(snp=rep(snps, 2),
                channel=rep(c("X", "Y"), each=length(snps)),
                as.data.frame(rbind(X, Y)))
rownames(result) <- 1:nrow(result)

# bring SNP rows together
result <- result[as.numeric(t(cbind(seq_along(snps), seq_along(snps)+length(snps)))),]
rownames(result) <- 1:nrow(result)

# write to fst file, maximally compressed
cat("Writing out FST...")
cat("Finished.\n")
dir.create("../output/objects")
write.fst(result, paste0("../output/objects/",prefix,"_intensities.fst"), compress=100)
rm(finrep)

#SNPs on sex X chromosome
Xsnp <- codes$marker[codes$chr == "X"]
Ysnp <- codes$marker[codes$chr == "Y"]

cat("\nGenerating Sex Intensity Diagnostic...\n")
#Pull X and Y Intensities using geneseek2qtl2.R method
#X intensities
Xres <- result[result$snp %in% Xsnp,]
Xint <- (Xres[Xres$channel == "X", -c(1,2)] + Xres[Xres$channel == "Y",-c(1,2)])/2
rownames(Xint) <- Xres$snp[Xres$channel == "X"]
#Y Intensities
Yres <- result[result$snp %in% Ysnp,]
rm(result)
Yint <- (Yres[Yres$channel == "X", -c(1,2)] + Yres[Yres$channel == "Y",-c(1,2)])/2
rownames(Yint) <- Yres$snp[Yres$channel == "Y"]

#T-test markers for significant intensity by sex
if(sexStat == "both"){
  #test X intensities
  xpval <- apply(t(Xint), 2, function(a) t.test(a~cross$covar[match(colnames(Xint),
                                                                    rownames(cross$covar)),
                                                              "sex"], 
                                                na.rm = T)$p.value)
  Xavg <- apply(Xint[xpval <= 0.05,], 2, mean)
  
  #test Y intensities
  ypval <- apply(t(Yint), 2, function(a) t.test(a~cross$covar[match(colnames(Yint),
                                                                  rownames(cross$covar)),
                                                            "sex"], 
                                              na.rm = T)$p.value)
  Yavg <- apply(Yint[ypval <= 0.05,], 2, mean)
}else{
  Xavg <- apply(Xint, 2, mean)
  Yavg <- apply(Yint, 2, mean)
}


intDF <- data.frame(Xavg = Xavg, Yavg = Yavg)

#Add in sex 
intDF$sex <- cross$covar[match(rownames(intDF),rownames(cross$covar)),"sex"]
intDF <- intDF[!is.na(intDF$sex),]
phetX <- rowSums(cross$geno$X == 2)/rowSums(cross$geno$X != 0)
intDF$phetX <- phetX[match(rownames(intDF),names(phetX))]

#Flag robust z-score outliers at double the suggested level of potential outliers
#Source: Boris Iglewicz and David Hoaglin (1993), 
#"Volume 16: How to Detect and Handle Outliers", 
#The ASQC Basic References in Quality Control: 
#Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.

flag <- character()
for(i in unique(intDF$sex)){
  tmp <- intDF[intDF$sex == i,]
  xthresh_hi <- median(tmp$Xavg) + (7*mad(tmp$Xavg))
  xthresh_lo <- median(tmp$Xavg) - (7*mad(tmp$Xavg))
  ythresh_hi <- median(tmp$Yavg) + (7*mad(tmp$Yavg))
  ythresh_lo <- median(tmp$Yavg) - (7*mad(tmp$Yavg))
  
  flag <- unique(c(flag, rownames(tmp[which(tmp$Xavg >= xthresh_hi |
                                              tmp$Xavg <= xthresh_lo |
                                              tmp$Yavg >= ythresh_hi |
                                              tmp$Yavg <= ythresh_lo),])))
  
}


a <- ggplot(intDF, aes(x = Xavg, y = Yavg, color = sex))+ 
  geom_point()+
  geom_label_repel(data = intDF[rownames(intDF) %in% flag,],
                   aes(label = rownames(intDF[rownames(intDF) %in% flag,])),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   show.legend = F)+
  labs(title = paste0("Sex Chromosome Diagnostic - ",prefix),
       x = "Average X Intensity", 
       y = "Average Y Intensity")+
  theme_bw()

#flagger for ggrepel here - swap out y for prop Het X
flag2 <- character()
for(i in unique(intDF$sex)){
  tmp <- intDF[intDF$sex == i,]
  xthresh_hi <- median(tmp$Xavg) + (7*mad(tmp$Xavg))
  xthresh_lo <- median(tmp$Xavg) - (7*mad(tmp$Xavg))
  phetXthresh_hi <- median(tmp$phetX) + (7*mad(tmp$phetX))
  phetXthresh_lo <- median(tmp$phetX) - (7*mad(tmp$phetX))
  
  flag2 <- unique(c(flag2, rownames(tmp[which((tmp$Xavg >= xthresh_hi |
                                                 tmp$Xavg <= xthresh_lo) | (
                                                   tmp$phetX >= phetXthresh_hi |
                                                     tmp$phetX <= phetXthresh_lo)),])))
  
}



b <- ggplot(intDF, aes(x = Xavg, y = phetX, color = sex))+ 
  geom_point()+
  labs(x = "Average X Intensity", 
       y = "Proportion of Heterogeneity on X Chr")+
  geom_label_repel(data = intDF[rownames(intDF) %in% flag2,],
                   aes(label = rownames(intDF[rownames(intDF) %in% flag2,])),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   show.legend = F)+
  theme_bw()



dir.create("../output/plots")
intDF <- saveRDS(intDF, paste0("../output/objects/",prefix,"_sexDiagData.rds"))

ggarrange(a,b, 
          ncol = 1, nrow = 2, 
          common.legend = T,
          legend = "right") %>% ggexport(filename = paste0("../output/plots/",prefix,"_SexDiag.pdf"))

sink("../output/sexDiag.txt")
cat("The following samples were flagged as outliers (robust Z score > 7)\nfrom the plot of Average X and Y Intensity:\n")
cat(flag, sep = "\n")
cat("The following samples were flagged as outliers (robust Z score > 7)\nfrom the plot of Average X Intensity and Proportion of Heterogeneity on chr X:\n")
cat(flag2, sep = "\n")
sink()

cat("fstGen.R: Status: Complete")