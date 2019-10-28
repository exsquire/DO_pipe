#Covariate check, visualizes potential batch effects in phenotypes data using covariate data 
covCheck <- function(cov, phen,
                     maxLev = 5,
                     pdfWid = 100, 
                     pdfHei = 100,
                     margins = c(100,80),
                     rowFont = 15,
                     colFont = 10,                  
                     path=paste0(getwd(),"/BatchEff_Pheno.pdf")){
  library(gplots)
  #cov should be a dataframe
  if(!is.data.frame(cov)){
    cov <- as.data.frame(cov)
  }
  
  #Ensure phen and covar samples are in the same order
  cov <- cov[match(rownames(phen), rownames(cov)),]
  
  if(all(rownames(cov) != rownames(phen))){
    stop("Warning! Row Mismatch Error.")
  }
  
  dumboDrop <- function(df){
    #If any variables are factors, convert them
    if(any(sapply(cov, is.factor))){
      cat("Converting factors to characters.\n")
      cov[,sapply(cov,is.factor)] <- sapply(cov[,sapply(cov,is.factor)],as.character)
    }

    #Transform covar into dummy variables 
    dumDown <- function(x, maxLev = 5){
      numLev <- length(unique(x))
      #Takes in a column, asks how many unique entries 
      if(numLev == 1 | numLev > maxLev){
        return(NULL)
      }else{
        dumDF <- data.frame(var = unique(x), dVar = seq(0,numLev - 1))
        
        #Forgive the weirdness, this gets around numeric covariates that would mess with dumDF 
        ind <- list()
        for(i in 1:numLev){
          #Build indexes that match dumDF key
          ind[[i]] <- which(x == dumDF[i,1])
        }
        for(j in 1:length(ind)){
          x[ind[[j]]] <- dumDF[j,2]
        }
        
        return(as.numeric(x))
      }
    }
    
    covarDum <- lapply(cov, dumDown, maxLev = maxLev)
    covarDum <- as.data.frame(covarDum[!sapply(covarDum,is.null)])
    
    return(covarDum)
  }
  
  dumCovar <- dumboDrop(cov)
  
  
  #Perform anova on each pheno for each covar
  batchCheck <- function(x, y){
    aovPval <- summary(aov(x ~ factor(y)))[[1]][1,5]
    return(aovPval)
  }
  
  outlist <- list()
  for(i in seq_along(colnames(dumCovar))){
    desig <- colnames(dumCovar)[i]
    tmp <- apply(phen, 2, batchCheck,y = dumCovar[,i])
    outlist[[desig]] <- tmp
  }
  
  #lapply pdjust
  outlist <- lapply(outlist, p.adjust, method = "BH")
  batchMat <- as.matrix(do.call("rbind", outlist))
  
  
  colBreaks = c(seq(0,0.05,length=1),
                seq(0.06,1,length=2))
  
  
  gradCol <- colorRampPalette(c("firebrick1",
                                "gray85"))(n=2)
  
  pdf(path,
      width = pdfWid, height = pdfHei)
  par(mar=c(1,1,1,1))
  heatmap.2(batchMat,Rowv=NA, Colv=NA, 
            dendrogram = "none",
            col = gradCol,
            breaks = colBreaks,
            scale = "none",
            #RowSideColors = RowSideColors,
            trace = "none",
            key=TRUE, 
            symkey = F,
            key.title = "",
            lmat = rbind(c(0,0,0,4,0,0),
                         c(0,1,1,1,1,1),
                         c(2,1,1,1,1,1),
                         c(0,1,1,1,1,1),
                         c(0,0,0,3,0,0)),
            lhei = c(0.2,0.5,0.5,0.5,0.1),
            lwid = c(0.1,2,2,2,2,2),
            margins = c(100,80),
            cexRow = rowFont,
            cexCol = colFont)
  dev.off()
}
