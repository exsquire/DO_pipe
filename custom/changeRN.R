#User controlled prefixing of rownames
changeRN <- function(x){
  if(any(is.null(rownames(x))) | any(is.na(rownames(x)))){
    stop("No rownames, try again.")
  }
  viewID <- head(rownames(x), n = 5)
  cat("First 5 Sample IDs...\n")
  cat(viewID, sep = "\n")
  altIDs <- "empty"
  while(!altIDs %in% c("Y","N")){
    altIDs <- toupper(
      readline("Would you like to alter sample IDs(Y/N)? "))
  }
  if(altIDs == "Y"){
    desig <- readline("Prefix sample IDs with...: ")
    samp <- paste0(desig, viewID)
    cat("IDs will now look like: \n")
    cat(samp, sep = "\n")
    conf1 <- toupper(readline("Confirm change (Y/N): "))
    
    if(conf1 == "Y"){
      rownames(x) <- paste0(desig, rownames(x))
      cat("IDs have been changed.\n")
    }else{
      cat("Aborting...\n
          Please try again.")
    }
    }else{
      cat("Returning unchanged.")
      return(x)
  }
  return(x)
}
