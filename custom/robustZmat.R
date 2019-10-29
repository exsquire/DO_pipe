#Use heatmap.2 to generate heatmap of robustZscores for each phenotype. Quickly spot outliers and their robust z-score positions across your phenotypes. 
#zoom = TRUE zooms in on rows and columns that contain robust Z scores > 7. 
robustZmat <- function(phen, 
                       prefix="",
                       rowFont = 2,
                       colFont =2,
                       margins = c(35,15),
                       pdfWid = 100,
                       pdfHei = 100,
                       path = getwd(),
                       zoom = FALSE,
                       valSize = 5){
  library(gplots)
  #Perform robust z-score transformation
  robZ <- function(x, y){
    zscoreMed <- (x[y]- median(x, 
                               na.rm = T))/mad(x, 
                                               na.rm = T)
    return(zscoreMed)
  }
  
  #Apply to get one subject
  outlist <- list()
  for(i in seq_along(rownames(phen))){
    desig <- rownames(phen)[i]
    tmp <- apply(phen, 2, robZ,y = i)
    outlist[[desig]] <- abs(tmp)
  }
  zMat <- as.matrix(do.call("rbind", outlist))
  
  #Clean for Inf values
  cat("Removing columns with Inf values.\n")
  zMat <- zMat[,!apply(zMat,2,function(x) any(is.infinite(x)))]
  if(zoom == FALSE){
    #We can set the RowSideColors using a vector, one initial fill followed by a which() or apply() statement to select the other color positions. 
    
    #https://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package - for formatting heatmaps
    
    #Add rowside colors for 0, < 5, <10, +10
    #RowSideColors <- rep("grey", nrow(zMat))
    #RowSideColors[apply(zMat > 7,1,any)] <- "firebrick2"
    
    #Take the range from 3.5 to max(zMat) and cut into 3 pieces
    outScale = quantile(seq(3.5, max(zMat, na.rm = T)))
    
    if(max(zMat, na.rm = T) < 3.5){
      stop("Robust z-score detects no outliers in data.")
    }
    
    colBreaks = c(seq(0,3.49,length=20),
                  seq(3.5,outScale[2],length=20),
                  seq(outScale[3],outScale[4],length=5),
                  seq(outScale[4]*1.01, max(zMat, na.rm = T),
                      length = 5))
    
    
    gradCol <- colorRampPalette(c("springgreen4", 
                                  "yellow1",
                                  "orange2",
                                  "firebrick1"))(n=length(colBreaks)-1)
    
    pdf(paste0(path,"/",prefix,"_RobZmat_Full.pdf"), width = pdfWid, height = pdfHei)
    par(mar=c(1,1,1,1))
    heatmap.2(zMat,Rowv=NA, Colv=NA, 
              dendrogram = "none",
              col = gradCol,
              breaks = colBreaks,
              #RowSideColors = RowSideColors,
              trace = "none",
              key=TRUE, 
              key.title = "",
              lmat = rbind(c(0,0,0,4,0,0),
                           c(0,1,1,1,1,1),
                           c(0,1,1,1,1,1),
                           c(0,1,1,1,1,1),
                           c(2,1,1,1,1,1),
                           c(0,1,1,1,1,1),
                           c(0,1,1,1,1,1),
                           c(0,0,0,3,0,0)),
              lhei = c(0.2,1,1,1,1,1,1,0.1),
              lwid = c(0.1,1,1,1,1,1),
              margins = margins,
              cexRow = rowFont,
              cexCol = colFont)
    dev.off()
  }else{
    if(max(zMat, na.rm = T) < 3.5){
      stop("Robust z-score detects no outliers in data.")
    }
    subZmat <- zMat[apply(zMat > 7,1,any, na.rm = T),
                    apply(zMat > 7,2,any, na.rm = T)]
    
    outScale = quantile(seq(3.5, max(subZmat, na.rm = T)))
    
    if(max(subZmat, na.rm = T) < 3.5){
      stop("Robust z-score detects no outliers in data.")
    }
    
    colBreaks = c(seq(0,3.49,length=20),
                  seq(3.5,outScale[2],length=20),
                  seq(outScale[3],outScale[4],length=5),
                  seq(outScale[4]*1.01, max(subZmat, na.rm = T),
                      length = 5))
    
    
    gradCol <- colorRampPalette(c("springgreen4", 
                                  "yellow1",
                                  "orange2",
                                  "firebrick1"))(n=length(colBreaks)-1)
    
    pdf(paste0(path,"/",prefix,"_RobZmat_Zoom.pdf"), width = pdfWid, height = pdfHei)
    par(mar=c(1,1,1,1))
    heatmap.2(subZmat,
              cellnote = round(subZmat,2), #add values
              notecol = "black", #value colors
              notecex = valSize, #value sizes
              Rowv=NA, Colv=NA, 
              dendrogram = "none",
              col = gradCol,
              breaks = colBreaks,
              #RowSideColors = RowSideColors,
              trace = "none",
              key=TRUE,
              key.title = "",
              lmat = rbind(c(0,0,0,4,0,0),
                           c(0,1,1,1,1,1),
                           c(0,1,1,1,1,1),
                           c(0,1,1,1,1,1),
                           c(2,1,1,1,1,1),
                           c(0,1,1,1,1,1),
                           c(0,1,1,1,1,1),
                           c(0,0,0,3,0,0)),
              lhei = c(0.2,1,1,1,1,1,1,0.1),
              lwid = c(0.1,1,1,1,1,1),
              margins = margins,
              cexRow = rowFont,
              cexCol = colFont)
    dev.off()
    
  }
}