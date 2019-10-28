library(qtl2)
library(qtlcharts)
library(broman)
library(kableExtra)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(gridExtra)

rm(list=ls())

#Recode function
allele_recode <- function(x){
  for(chr in seq_along(x$founder_geno)) {
    fg <- x$founder_geno[[chr]]
    g <- x$geno[[chr]]
    f1 <- colSums(fg==1)/colSums(fg != 0)
    
    
    fg[fg==0] <- NA
    g[g==0] <- NA
    
    fg[,f1 < 0.5] <- 4 - fg[,f1 < 0.5]
    g[,f1 < 0.5]  <- 4 - g[,f1 < 0.5]
    
    fg[is.na(fg)] <- 0
    g[is.na(g)] <- 0
    
    x$founder_geno[[chr]] <- fg
    x$geno[[chr]] <- g
  }
  return(x)
}
#-------------------------------------------------------
cat("\nCalling script: proc_cross.R\n")

#Run script from folder with the 3 input files 
#Form a list of filepaths where [[1]] == cross2, 
#[[2]] == codes, and [[3]] == Final Report
pathIn <- list(
  crossPath = list.files("../../inputs/", pattern = "cross", full.names = T, ignore.case = T),
  prefPath = list.files("../../inputs/", pattern = "prefix", full.names = T, ignore.case = T)
)

#Check
if(!all(sapply(pathIn, function(x) length(x) > 0))){
  cat("Error: Input could not be found.\n")
  print(sapply(pathIn, function(x) length(x) > 0))
}else{
  cat("\nFound the following inputs: \n")
  print(sapply(pathIn, function(x) length(x) > 0))
}

#Inputs 
input <- input0 <- readRDS(pathIn[[1]])
prefix <- readRDS(pathIn[[2]])

#Outputs
msgOutputs <- list()

#Quick Clean-up
##Drop null markers
cat("\nDropping null markers...\n")
input <- drop_nullmarkers(input)
#drop null output
null0 <- sum(sapply(input0$geno,ncol, simplify = T)) - 
  sum(sapply(input$geno,ncol, simplify = T))
rm(input0)
msgOutputs[["NullMarks"]] <- paste0("Dropped ", null0, " markers.")



##Recode Genotypes - swap major and minor allele factor when there are more minors than majors
cat("\nRecoding alleles to reflect 'Major' and 'Minor' Status\n")
input <- allele_recode(input)
#-----------------------------------------

cat("\nGenerating Missing Genotype Diagnostic files...\n")

percent_missing <- n_missing(input, "ind", "prop")*100

labels <- paste0(names(percent_missing), " (", round(percent_missing), "%)")


percMissDF <- data.frame(samp = labels,
                         Index = seq_along(labels), 
                         val = round(percent_missing,2))
saveRDS(percMissDF, file = paste0("../output/objects/",prefix,"_percMiss_SampDF.rds"))

pc_thresh <- median(percMissDF$val)+(3.5*mad(percMissDF$val))

#GGPLOT
percMiss <- ggplot(percMissDF, 
                   aes(x =Index, y = val))+
  geom_point(color = "dodgerblue3", 
             size = 3)+
  scale_y_continuous(limits = c(0, 100))+
  labs(title = "Missing Genotype Diagnostic",
       y = "Missing Genotype Data (%)",
       x = "Mouse")+ 
  geom_label_repel(data = percMissDF[percMissDF$val >= pc_thresh,],
                   aes(label = samp),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   ylim = c(max(percMissDF$val) + 5, NA),
                   segment.color = 'grey50')+
  theme_bw()

percMiss %>% ggexport(filename = paste0("../output/plots/",prefix,"_PercMiss_Plot.pdf")) 
###################################
#Check for duplicate Mice - output suspicious ones and rug hist
cat("\nChecking for duplicate samples...\n")
cat("\nRunning pairwise genotype scan...")
cg <- compare_geno(input, cores=10)
cat("Finished.\n")

cat("\nGenerating Duplicate Sample Diagnostic...")
summary(cg)
summary(cg[upper.tri(cg)])
msgOutputs[["DupesMsg"]] <- summary(cg)
msgOutputs[["Dupes_Summary"]]<- summary(cg[upper.tri(cg)])

#No duplicates found - in the plot below, we see no ticks above the 0.9 mark. 
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
                     base_size = 15,
                     padding = unit(c(2, 4), "mm"))
tbl <- tableGrob(matrix(round(summary(cg[upper.tri(cg)]),3),nrow = 1), 
                 rows=NULL, 
                 cols = c("Min","1st Qu.","Median","Mean","3rd Qu.","Max"), 
                 theme=tt)

dupePlot <- ggplot(data.frame(props = cg[upper.tri(cg)]), aes(x = props)) + 
  geom_histogram(breaks=seq(0, 1, length=201)) + 
  scale_x_continuous(limits = c(0,1)) + 
  theme_bw()+
  geom_rug(colour = "dodgerblue3") + 
  geom_vline(xintercept =  0.9,
             colour = "firebrick2",
             linetype = "dashed" ) + 
  labs(title = "Duplicate Sample Diagnostic",
       x= "Proportion of Matching Genotypes",
       y= "Counts")

grid.arrange(dupePlot, tbl, nrow = 2, as.table = T, heights = c(2,1)) %>%
  ggsave(filename = paste0("../output/plots/",prefix,"_DupeDiag.pdf"), width = 6, height = 4)

#If large number of samples with missing data are present, remove to denoise
#cgsub <- cg[percent_missing < 50, percent_missing < 50]
#par(mar=c(5.1,0.6,0.6, 0.6))
#hist(cgsub[upper.tri(cgsub)], breaks=seq(0, 1, length=201),
#     main="", yaxt="n", ylab="", xlab="Proportion matching genotypes")
#rug(cgsub[upper.tri(cgsub)])
cat("Finished.\n")

#-------------------------------------------------------------
#Genotype Tri-plot
g <- do.call("cbind", input$geno[1:19])
fg <- do.call("cbind", input$founder_geno[1:19])

#Subset for columns that do not have any zero genotypes. 
g <- g[,colSums(fg==0)==0]
fg <- fg[,colSums(fg==0)==0]
cat("\nExporting 'g' and 'fgn' objects...")
saveRDS(g, paste0("../output/objects/",prefix,"_g.rds"))
#Vector for number of minor alleles by SNP
fgn <- colSums(fg==3)
saveRDS(fgn, paste0("../output/objects/",prefix,"_fgn.rds"))
cat("Finished.\n")

#Calc genofreqs by individual
gf_ind <- vector("list", 4)
for(i in 1:4) {
  gf_ind[[i]] <- t(apply(g[,fgn==i], 1, function(a) table(factor(a, 1:3))/sum(a != 0)))
}
cat("\nGenerating Sample Tri-plot...")
#Broman: The following triangle plots show the genotype frequency distributions for the mice, among the four groups of markers with common minor allele frequency (MAF) in the founder strains. These plots make use of the fact that for a point within an equilateral triangle, the sum of the distances to the three sides is a constant.
pdf(paste0("../output/plots/",prefix,"_SampTriPlot.pdf"), width = 6, height = 6)
par(mfrow=c(2,2), mar=c(0.6, 0.6, 2.6, 0.6))
for(i in 1:4) {
  triplot(c("AA", "AB", "BB"), main=paste0("MAF = ", i, "/8"))
  tripoints(gf_ind[[i]], pch=21, bg="lightblue")
  tripoints(c((1-i/8)^2, 2*i/8*(1-i/8), (i/8)^2), pch=21, bg="violetred")
  
  if(i>=3) { # label mouse with lowest het
    wh <- which(gf_ind[[i]][,2] == min(gf_ind[[i]][,2]))
    tritext(gf_ind[[i]][wh,,drop=FALSE] + c(0.02, -0.02, 0),
            names(wh), adj=c(0, 1))
  }
  
  # label other mice
  if(i==1) {
    lab <- rownames(gf_ind[[i]])[gf_ind[[i]][,2]>0.3]
  }
  else if(i==2) {
    lab <- rownames(gf_ind[[i]])[gf_ind[[i]][,2]>0.48]
  }
  else if(i==3) {
    lab <- rownames(gf_ind[[i]])[gf_ind[[i]][,2]>0.51]
  }
  else if(i==4) {
    lab <- rownames(gf_ind[[i]])[gf_ind[[i]][,2]>0.6]
  }
  
  for(ind in lab) {
    if(grepl("^F", ind) && i != 3) {
      tritext(gf_ind[[i]][ind,,drop=FALSE] + c(-0.01, 0, +0.01), ind, adj=c(1,0.5))
    } else {
      tritext(gf_ind[[i]][ind,,drop=FALSE] + c(0.01, 0, -0.01), ind, adj=c(0,0.5))
    }
  }
}
dev.off()
cat("Finished.\n")
#-------------------------------
#Generate Cross_Pass1
#-------------------------------
cat("\nFiltering Cross Raw for missing samples...\n")

#At this point, kick out a cross_pass1, which 
#removes samples with large +3.5mad from median missing or > 20%
#And dropped all null markers
pc_thresh2 <- median(percMissDF$val)+(3.5*mad(percMissDF$val))
if(pc_thresh2 < 20){
  cross_pass1 <- input[percent_missing < pc_thresh2]
}else{
  cross_pass1 <- input[percent_missing < 20]
}
cat("Missing percent threshold:", round(pc_thresh2,2))
cat("\nExporting Cross Pass1 object...")
saveRDS(cross_pass1, file = paste0("../output/objects/",prefix,"_cross_pass1.rds"))
cat("Finished.\n")
#--------------------------------
#Generate missing data by Marker plot - possible after pass1
cat("\nGenerating Global Marker Diagnostic...")
pmis_mar <- n_missing(cross_pass1, "marker", "proportion")*100

tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
                     base_size = 10,
                     padding = unit(c(2, 4), "mm"))
tbl <- tableGrob(matrix(round(summary(pmis_mar),3),nrow = 1), 
                 rows=NULL, 
                 cols = c("Min","1st Qu.","Median","Mean","3rd Qu.","Max"), 
                 theme=tt)

cuts <- numeric()
for(i in seq(0,90,by = 10)){
  cuts <- c(cuts, round(sum(which(pmis_mar > i))/sum(which(pmis_mar > 0))*100,2))
}

tbl2 <- tableGrob(matrix(cuts,nrow = 1), 
                 rows=NULL, 
                 cols = paste0(">",seq(0,90,10),"%"), 
                 theme=tt)

misMarPlot <- ggplot(data.frame(props = pmis_mar), aes(x = props)) + 
  geom_histogram(breaks=seq(0, 100, length=201)) + 
  scale_x_continuous(limits = c(0,100)) + 
  theme_bw()+
  geom_rug(colour = "dodgerblue3") + 
  labs(title = "Missing Genotypes - Marker Diagnostic",
       x= "Missing Genotypes (%)",
       y= paste0("Counts (n = ",length(pmis_mar),")"))

grid.arrange(misMarPlot, tbl, tbl2,nrow = 3, as.table = T, heights = c(4,1,1)) %>%
  ggsave(filename = paste0("../output/plots/",prefix,"_MissMarkDiag.pdf"), width = 6, height = 4)
cat("Finished.\n")

invisible(gc())
cat("\nproc-cross.R: Status: Complete")























