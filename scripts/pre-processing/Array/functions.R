# filtering functions
filt.by.neg.probes <- function(data, cdata, qdist = .9, perc = 50, verbose = TRUE) {
  if(verbose) cat(" - filtering by control-negative probes..\n")
  if(verbose) print(paste0("Qdist - ", qdist, ", Perc - ", perc))
  # compile a threshold value from the negative control probes for each sample
  thres.vec <- apply(cdata, 2, function(x) quantile(x, qdist, na.rm=TRUE))
  # filter out the probes whose...
  score.neg <- matrix(0, nrow=nrow(data), ncol=ncol(data))
  for(i in 1:dim(data)[2])  score.neg[which(data[,i] > thres.vec[i]),i] <- 1
  neg.sumrow <- apply(score.neg, 1, sum)
  if(verbose) print(paste0("Rounded - ", round(ncol(data)*perc/100)))
  filt.data <- data[neg.sumrow >= round(ncol(data)*perc/100),]
  if(verbose) cat("     * number of probes to remove:", (nrow(nc.data) - nrow(filt.data)), "\n")
  return(filt.data)
}

filt.by.var.dup.probes <- function(data, qdist = .9, perc = 50, verbose = TRUE, plot = TRUE, ...) {
  if(verbose) cat(" - filtering by variance duplicated probes..\n")
  # build the list of indexes for each duplicated non control probe
  un.probes <- unique(rownames(data))
  list.id.dup.probes <- lapply(un.probes, function(p) {which(rownames(data) == p)})
  un.probes <- un.probes[-which(sapply(list.id.dup.probes, length) < 2)]
  list.id.dup.probes <- list.id.dup.probes[-which(sapply(list.id.dup.probes, length) < 2)]
  score.neg <- matrix(0, nrow=length(list.id.dup.probes), ncol=ncol(data))
  for(j in 1:ncol(data)) {
    dist.var <- unlist(lapply(list.id.dup.probes, function(i) var(data[i,j])))
    thres.vec <- quantile(dist.var, qdist, na.rm=TRUE)
    score.neg[which(dist.var > thres.vec), j] <- 1
  }
  neg.sumrow <- apply(score.neg, 1, sum)
  probs.to.rem <- un.probes[neg.sumrow >= round(ncol(data)*perc/100)]
  probs.to.rem <- which((rownames(data) %in% probs.to.rem))
  if(verbose) cat("     * number of probes to remove:", length(probs.to.rem), "\n")
  return(data[-probs.to.rem,])
}

filtering_negative_duplicate = function(grid, matrix_intensity, nc.data){
  filt.data <- apply(grid, 1, FUN = function(x) {
    if(verbose) cat("*** (qidst:",x[1],", perc.samples:",x[2], ")\n", sep="")
    fd <- filt.by.neg.probes(data = matrix_intensity, nc.data, qdist = x[1], perc = x[2])
    filt.by.var.dup.probes(fd, qdist = x[1], perc = x[2], verbose=TRUE)
  })
  return(filt.data)
}

select_platform <- function(dataset_ensembl, platform){
  # platform = 'agilent' or 'affy'
  ensembl <- useMart('ensembl', dataset = dataset_ensembl)
  tables <- listAttributes(ensembl)
  
  if(nrow(tables) > 0){
    print(tables[grep(platform, tables[,1]),]) 
  }else{
    stop("Platform not available, however if you meant affymetric use 'affy' and 'agilent' for agilent related platforms")
  }
}

get_annotation <- function(platform, probes, dataset_ensembl){
  
  ensembl <- useMart('ensembl', dataset = dataset_ensembl)
  annot <- getBM(
    attributes = c(platform,
                   'wikigene_description',
                   'ensembl_gene_id',
                   'entrezgene_id',
                   'external_gene_name'),
    filters = platform,
    values = probes,
    mart = ensembl)
  
  annot <- merge(
    x = as.data.frame(probes),
    y =  annot,
    by.y = platform,
    all.x = T,
    by.x = 'probes')
  
  return(annot)
}

harmonise <- function(pd){
  # harmonising metadata table
  if(!"group" %in% colnames(pd)){
    colnames(pd)[colnames(pd) == "exposure"] = "group"
  }
  
  pd$title = as.character(as.vector(pd$title))
  
  if(any(grepl("Sample FAKE", pd$title))){
    pd = pd[!pd$title == "Sample FAKE",]
  }
  
  # change label_ch1 to dye
  colnames(pd)[colnames(pd) == "label_ch1"] = "dye"
  return(pd)
}

load_phenotype <- function(inputfile) {
  tryCatch({
    
    phTable <- readxl::read_xlsx(inputfile, col_names = TRUE)
    
    # save original table
    orig_phTable <- phTable
    
    # cleaning and conversions
    coltypes <- sapply(phTable, class)
    coltypes.logical.idx <- which(coltypes == "logical")
    if (length(coltypes.logical.idx) > 0) {
      phTable[coltypes.logical.idx] <- lapply(phTable[coltypes.logical.idx], as.character)
    }
    
    # remove columns with only 1 level
    nrlevels <- sapply(phTable, function(x) length(unique(x)))
    phTable <- phTable[, nrlevels > 1, drop = FALSE]
    
    # creation of the class table 
    phClassTable <- data.frame(
      Variable = names(phTable),
      Class = sapply(phTable, class),
      Type = ifelse(sapply(phTable, is.character), "factor", "vector")
    )
    
    # return the results 
    return(list(orig_phTable, phTable, phClassTable))
    
  }, error = function(e) {
    message("Error: ", e$message)
    NULL
  })
}


convert_to_chr <- function(phClassTable){
  factorIdx <- which(phClassTable$Type=="factor")
  factorCols <- NULL
  if(length(factorIdx)>0){
    factorCols <- colnames(phTable)[factorIdx]
    phFactor <- factorize_cols(phTable, factorIdx)
  }else{
    phFactor <- phTable
  } 
  return(phFactor)
}

upload_rawdata <- function(phTable, orig_phTable,  phFactor, removedSamplesInfo,
                              fileNameColName, sampleColName, dyeColName,
                              celDir, arrType, affCDF
                              ){


    phTable <- clean_phTable
    orig_phTable <- orig_phTable
    phFactor <- phFactor

    removedSamplesInfo <- removedSamplesInfo
    fileNameColName <- fileNameColName
    sampleColName <- sampleColName
    dyeColName <- dyeColName
    fileNameColID <- which(colnames(phTable) %in% fileNameColName)
    sampleColID <- which(colnames(phTable) %in% sampleColName)
    dyeColID <- which(colnames(phTable) %in% dyeColName)

    celDir <- celDir
    arrType <- arrType
    affCDF <- affCDF

    #fileNames <- unique(phTable[,fileNameColName])    
    fileNames <- unique(as.data.frame(orig_phTable)[,fileNameColName]) %>%
      sub(".*/([^/]+)$", "\\1", .)#sub(".*/([^/]+)$", "\\1")
    fileNameCount <- length(fileNames)

    #Check missing files
    if(arrType=="il_methyl"){    #riflessioni123    
      # fileNamesTmp <- as.vector(sapply(fileNames, function(x){g<-paste0(x, "_Grn.idat"); r<-paste0(x, "_Red.idat"); return(c(g,r))}))
      # fileNameCount <- length(fileNamesTmp)
      # fileChk <- which(fileNamesTmp %in% dir(celDir))
      # fileChkCount <- length(fileChk)
      # missingFiles <- fileNamesTmp[-fileChk]
      # missingFilesStr <- paste0(missingFiles, collapse="\n")
    }else{
      fileChk <- which(fileNames %in% dir(celDir))
      fileChkCount <- length(fileChk)
      if(length(fileChk)>0){
        missingFiles <- fileNames[-fileChk]
      }else{
        missingFiles <- fileNames
      }
      missingFilesStr <- paste0(missingFiles, collapse="\n")
    }

    if(fileChkCount < fileNameCount){

        stop(paste0("Could not find the needed raw data files in the selected directory!\n\nPlease check and select the correct directory containing the raw data files specified in the phenotype data.\n\nFound: ", fileChkCount, " of ", fileNameCount, ". \n\nMISSING FILES:\n\n", missingFilesStr))
    }

    if(arrType=="af_exp"){
      #checks presence of gz file to unzip
      fileNames <- unique(as.data.frame(orig_phTable)[,fileNameColName]) %>%
        sub(".*/([^/]+)$", "\\1", .)#sub(".*/([^/]+)$", "\\1")

      tmp_extensions <- tools::file_ext(orig_phTable$supplementary_file) %>% unique %>% tolower(.)

      if (any (tmp_extensions == "gz") ){
        phTable$supplementary_file <-    sapply( seq_along(fileNames), function(i)  {
          if(tools::file_ext(fileNames[i]) %in% "gz") {
            tryCatch({
                tmp_celfilename <- file.path(celDir, fileNames[i]) %>%
                                      R.utils::gunzip(., overwrite = TRUE, remove = FALSE)

               x<- as.character(tmp_celfilename[1]) %>%
                  sub(".*/([^/]+)$", "\\1", .)
               x
            },
            error = function(e) {
              message("Error unzipping file at ", i, "th position - ", e$message)
            }
            )
          }
        })
      }

      fileNames <- unique(phTable$supplementary_file)
      tmpFile <- fileNames[which(tools::file_ext(fileNames) %in% c("cel", "CEL"))[1]] %>%
        file.path(celDir, .)


      #Get the CDF name from a user provided CEL file
      celHeader <- affyio::read.celfile.header(tmpFile)

      cdfName <- gsub("-|_", "", celHeader$cdfName)

      #Get name of all installed R packages
      allInstalledPkgs <- rownames(installed.packages())

      #Format the CDF name to Bioconductor format
      cdfNameBioc <- paste0(tolower(cdfName), "cdf")
      print("Bioconductor cdfName:")
      print(cdfNameBioc)

      qcCDF <- affCDF
      warned <- FALSE

      #Check the bioconductor library for CDF is installed. If yes then call library() else try to install the associated Bioconductor package
      if(!any(grepl(cdfNameBioc, allInstalledPkgs))){
        tryCatch(BiocManager::install(cdfNameBioc)) #BiocInstaller::biocLite(cdfNameBioc, suppressUpdates=TRUE))
        allInstalledPkgs <- rownames(installed.packages())
        if(any(grepl(cdfNameBioc, allInstalledPkgs))){
          qcCDF <- cdfNameBioc
          library(package=qcCDF, character.only=TRUE)
        }
      }else{
        qcCDF <- cdfNameBioc
        library(package=qcCDF, character.only=TRUE)
      }

      #Check if user selected CDF annotation mathces with the CDF name in the affymetrix CEL file
      if(!is.null(affCDF) & !any(grepl(cdfName, affCDF, ignore.case=T))){
       stop(
       paste0("CDF annotation mismatch with the CEL files!\n\nNeed CDF file corresponding to '", cdfName, "'\n\nUser provided CDF '", affCDF, "'")
         )
      }

      rownames(phTable) <- phTable[,sampleColName]
      rownames(phFactor) <- phTable[,sampleColName]
      pheno <- new("AnnotatedDataFrame", data=phFactor)

      #Read affymetrix cel files
      eset <- affy::justRMA(filenames=fileNames, celfile.path=celDir, phenoData=pheno, sampleNames=phTable[,sampleColName], normalize=TRUE, background=TRUE, cdfname=affCDF)
      exprs <- exprs(eset)

      #Check if QC can be performed for the provided annotation CDF
      print("qcCDF:")
      print(qcCDF)
      #Create affyBatch object for QC
      err <- 0
      tryCatch(setQCEnvironment(qcCDF), error=function(e){err<<-1})
      if(!err){
        celFilePaths <- file.path(celDir, fileNames)
        affyBatchObject <- affy::read.affybatch(filenames=celFilePaths, phenoData=pheno, cdfname=qcCDF)
      }



      #return multiple outputs
      return (list (phTable = phTable, affyBatchObject = affyBatchObject, exprs = exprs,
                    qcCDF = qcCDF, err= err ) )



    }else if(arrType=="il_methyl"){
      pdFactor <- phFactor
      pdFactor$Basename <- factor(file.path(celDir, fileNames))
      RGset <- minfi::read.metharray.exp(targets=pdFactor)
      detP <- minfi::detectionP(RGset)

      #Write output
      RGsetFile <- file.path(dataDir, "RGset.rds")
      saveRDS(object=RGset, file=RGsetFile)
      taskManager$addOutput('phTableFile', 'file', phTableFile)

      taskManager$addOutput('detP', 'numeric', detP)
    }
    loadedRaw <- TRUE

    taskManager$addOutput('loadedRaw', 'logical', loadedRaw)
    ### END ###

    # Add an error
    taskManager$addError('error1', "Error during x calculation")
    taskManager$end()



}

unzip_gz_files <- function(file_column, filedir) {
  gz_files <- file_column[grepl("\\.gz$", file_column)]  %>% file.path(filedir, .)
  
  for (file in gz_files) {
    message("Unzipping: ", file)
    tryCatch({
      R.utils::gunzip(file, overwrite = TRUE, remove = FALSE)
    }, error = function(e) {
      message("Error unzipping file: ", file, " - ", e$message)
    })
  }
}

## PRE-PROCESSING

get.info <- function(pd, verbose = TRUE, ...) {
  vars <- colnames(pd)
  print(data.frame(Variables = vars))
  n1 <- NULL
  while(is.null(n1)) {  
    n1 <- read.input.number(" - indicate the variable of interest:")
    if(!is.numeric(n1)) n1 <- NULL
    if(n1 > length(vars)) n1 <- NULL
  }
  var.int <- vars[n1]
  vars <- vars[-n1]
  # n1 <- read.input.number(" - how many contrasts:")
  levs <- names(table(pd[,var.int]))
  print(data.frame(Levels = levs))
  contrasts <- c("")
  flag <- TRUE
  while(flag) {
    cont <- read.input.word(message = " - indicate a contrast:")
    cont <- as.numeric(unlist(strsplit(cont, ",")))
    contrasts <- c(contrasts, paste(levs[cont[1]],levs[cont[2]],sep="-")) 
    flag <- continue.read.inputs()
  }
  contrasts <- contrasts[-1]
  n2 <- NULL
  # indicating the covariates
  covariates <- NULL
  while(!is.numeric(n2)) {
    print(data.frame(Variables = vars))
    n2 <- read.input.number(" - add a covariate <enter '0' to go to the next step>:")
    if(n2 == 0) {
      n2 <- NULL
      break
    }
    if(is.null(covariates))  covariates <- n2
    else  covariates <- c(covariates, n2)
    n2 <- NULL
  }
  if(!is.null(covariates)) {
    covs <- vars[covariates]
    vars <- vars[-covariates]
  }
  else  covs <- NULL
  # indicating the batch variables 
  batches <- NULL
  while(!is.numeric(n2)) {
    print(data.frame(Variables = vars))
    n2 <- read.input.number(" - add a batch <enter '0' to go to the next step>:")
    if(n2 == 0) break
    if(is.null(batches)) batches <- n2
    else batches <- c(batches, n2)
    n2 <- NULL
    print(batches)
  }
  if(!is.null(batches)) bats <- vars[batches]
  else bats <- NULL
  list(var.int = var.int,
       covariates = covs,
       batches = bats,
       contrasts = contrasts)
}

read.contrasts <- function(group) {
  levs <- names(table(group))
  print(data.frame(Levels = levs))
  contrasts <- c("")
  flag <- TRUE
  while(flag) {
    cont <- read.input.word(message = " - indicate a contrast:")
    cont <- as.numeric(unlist(strsplit(cont, ",")))
    contrasts <- c(contrasts, paste(levs[cont[1]],levs[cont[2]],sep="-")) 
    flag <- continue.read.inputs()
  }
  contrasts <- contrasts[-1]
  contrasts
}

run.norm <- function(filt.data, rgList, method, method2, plot=FALSE){
  norm.data <- NULL
  
  if(plot)  {
    boxplot(log2(filt.data), las=2, cex=0.7, main="Before normalization.")
  }
  
  if(is.matrix(filt.data) || is.data.frame(filt.data)){
    print("Data is in matrix format. Performing log2 transformation before normalization")
    filt.data <- log2(filt.data)
  }
  
  norm_switch <- function(method){
    switch(method,
           "quantile" = {
             if(is.matrix(filt.data)){
               limma::normalizeQuantiles(filt.data)
             }else{
               stop("Input is not a numeric matrix!")
             }
           },
           "vsn" = limma::normalizeVSN(filt.data),
           "cl" = {
             if(is.matrix(filt.data) || is.matrix(as.matrix(filt.data))){
               limma::normalizeCyclicLoess(filt.data, method="fast")
             }else{
               stop("Input is not a numeric matrix nor could be coerced to a matrix!")
             }
           },
           "BA" = limma::normalizeBetweenArrays(filt.data, method=method2)
    )
  }
  #norm.data <- normalizeQuantiles(log2(filt.data)) 
  norm.data <- norm_switch(method) 
  
  if(plot)  {
    boxplot(norm.data, las=2, cex=0.7, main="After normalization.")
  }
  
  return(norm.data)
} 

pre.proc <- function(pd, rgList, var.mds, qdist = c(.75, .9), perc = c(75, 50), verbose = TRUE, plot = TRUE, ...) {
  norm.data <- NULL
  data <- cbind(rgList$G, rgList$R)
  rownames(data) <- rgList$genes$ProbeName
  # defining the colnames for data
  print(head(colnames(data)))
  print(head(rownames(pd)))
  colnames(data) <- rownames(pd)
  nc.data <- data[which(rgList$genes$ControlType==0),]
  c.data  <- data[which(rgList$genes$ControlType==-1),]
  # filter with different combinations of qdist and perc
  if(length(qdist) != 1 & length(perc) != 1) {
    grid <- expand.grid(qdist, perc)
    filt.data <- apply(grid, 1, FUN = function(x) {
      if(verbose) cat("*** (qidst:",x[1],", perc.samples:",x[2], ")\n", sep="")
      #fd <- 
      filt.by.neg.probes(nc.data, c.data, qdist = x[1], perc = x[2], verbose)
      #filt.by.var.dup.probes(fd, qdist = x[1], perc = x[2], verbose)
    })
    # print out some information
    grid.info <- apply(grid, 1, function(x) paste(x[1], x[2], sep="-"))
    perc.probes <- unlist(lapply(filt.data, function(f) {
      round(100 - (((nrow(nc.data) - nrow(f))/nrow(nc.data))*100), digits = 2)}))
    print(data.frame(grid.info, perc.probes, num.probes = unlist(lapply(filt.data, nrow))))
    #n1 <- read.input.number(message = paste("Insert a number between 1 and", length(grid.info)), more.info=NULL)
    n1 <- 1 
    # normalize the data
    if(verbose) cat(" - quantile normalization method..\n")
    if(plot)  {
      boxplot(log2(filt.data[[n1]]), las=2, cex=0.7, main="Before normalization.")
    }
    norm.data <- normalizeQuantiles(log2(filt.data[[n1]]))
    if(plot)  {
      boxplot(norm.data, las=2, cex=0.7, main="After normalization.")
    }
  }
  else {
    if(verbose) cat("*** (qdist:",qdist,", perc.samples:",perc, ")\n", sep="")
    filt.data <- filt.by.neg.probes(nc.data, c.data, qdist = qdist, perc = perc, verbose)
    #filt.data <- filt.by.var.dup.probes(filt.data, qdist = qdist, perc = perc, verbose)
    cat("perc.probes", round(100 - (((nrow(nc.data) - nrow(filt.data))/nrow(nc.data))*100), digits = 2), "\n")
    cat("num.probes", nrow(filt.data), "\n")
    if(plot)  {
      boxplot(log2(filt.data), las=2, cex=0.7, main="Before normalization.")
    }
    norm.data <- normalizeQuantiles(log2(filt.data)) 
    if(plot)  {
      boxplot(norm.data, las=2, cex=0.7, main="After normalization.")
    }
  }
  return(norm.data)
} 



eval.by.neg.probes <- function(data, cdata, qdist = .9, perc = 50, verbose = TRUE, ...) {
  if(verbose) cat(" - evaluating by control-negative probes..\n")
  if(verbose) print(paste0("Qdist - ", qdist, ", Perc - ", perc))
  # compile a threshold value from the negative control probes for each sample
  thres.vec <- apply(cdata, 2, function(x) quantile(x, qdist, na.rm=TRUE))
  # filter out the probes whose...
  score.neg <- matrix(0, nrow=nrow(data), ncol=ncol(data))
  for(i in 1:dim(data)[2])  score.neg[which(data[,i] > thres.vec[i]),i] <- 1
  neg.sumrow <- apply(score.neg, 1, sum)
  if(verbose) print(paste0("Rounded - ", round(ncol(data)*perc/100)))
  eval.vector <- neg.sumrow >= round(ncol(data)*perc/100)
  if(verbose) cat("     * number of poor quality probes:", (nrow(data) - sum(eval.vector)), "\n")
  return(eval.vector)
}



## VISUALIZE/IDENTIFY & REMOVE BATCH EFFECTS
monitor.technical.variation <- function(data, pd, npc = 10, verbose = TRUE, ...) {
  test <- sapply(colnames(pd), function(b) length(table(pd[,b])) > 1 & length(table(pd[,b])) != length(pd[,b]))
  # generate the counfounding plot
  print(test)
  res <- confounding(pd[,test], margins = c(10,10))
  pr <- prince(data, pd[,test], top=npc)
  print(pr)
  # generate the prince plot
  prince.plot(prince=pr, margins = c(15,15))
  # hc plot
  hca.plot(data, pd[,test], method = "correlation")
  list(conf = res, pr = pr)
}

pc.anaylsis <- function(data, pd, npc = 10, verbose = TRUE, ...) {
  print("PC-analysis")
  # run the principal-component analysis
  var.names <- colnames(pd)
  # remove variables with one level
  test <- sapply(var.names, function(v) { 
    if(is.factor(pd[,v]))
      if(nlevels(pd[,v]) == 1) TRUE
    else FALSE
    else FALSE}) 
  ids <- which(test == TRUE)
  if(length(ids) > 0) var.names <- var.names[-ids]
  # run princomp
  pc <- princomp(data, cor=T)
  eig <- with(pc, unclass(loadings))
  if(ncol(eig) > npc) eig <- eig[,1:npc]
  perc.importance <- round(pc$sdev^2/sum(pc$sdev^2)*100, 2)
  # compile assoc matrix
  assoc <- matrix(NA, nrow=length(var.names), ncol=ncol(eig))
  rownames(assoc) <- var.names
  colnames(assoc) <- colnames(eig)
  for(i in 1:npc) {
    assoc[,i] <- sapply(var.names, function(a) {
      if(is.factor(pd[,a])) 
        anova(lm(eig[,i]~as.factor(pd[,a])))[1,5]
      else if(is.numeric(pd[,a]))
        anova(lm(eig[,i]~pd[,a]))[1,5]
    })
  }
  colnames(assoc) <- paste(rep("PC", npc), 1:npc, " (", perc.importance[1:npc], "%) ", sep="")
  return(assoc)
}

pc.anaylsis.2 <- function (g, o, npc = 10, center = T){
  test <- sapply(colnames(o), function(b) length(table(o[,b])) > 1 & length(table(o[,b])) != length(o[,b]))
  if(nrow(o) == 1) {
    o <- as.data.frame(o)
    colnames(o) <- "sva.1"
  }
  print(test)
  if (is.matrix(g) != T) {
    stop("g is not a matrix")
  }
  if (is.data.frame(o) != T) {
    stop("o is not a data.frame")
  }
  classes <- unlist(lapply(unclass(o), class))
  if (any(classes == "character")) {
    stop("o contains characters")
  }
  nrlevels <- unlist(lapply(unclass(o), function(x) length(levels(x))))
  if (any(nrlevels == 1)) {
    stop("o contains factors with only one level")
  }
  if (npc > ncol(g)) {
    stop("npc is larger than ncol(g)")
  }
  if (npc > nrow(g)) {
    stop("npc is larger than nrow(g)")
  }
  if (!identical(rownames(o), colnames(g))) {
    warning("Colnames of g are not the same as rownames of o")
  }
  if (center == T) {
    pr <- prcomp(t(g))
  }
  if (center == F) {
    pr <- prcomp(t(g), center = F)
  }
  linp <- matrix(ncol = npc, nrow = ncol(o))
  rownames(linp) <- colnames(o)
  rsquared <- matrix(ncol = npc, nrow = ncol(o))
  rownames(rsquared) <- colnames(o)
  for (i in 1:ncol(o)) {
    for (j in 1:npc) {
      fit <- lm(pr$x[, j] ~ o[, i])
      s <- summary(fit)
      linp[i, j] <- pf(s$fstatistic[1], s$fstatistic[2], 
                       s$fstatistic[3], lower.tail = FALSE)
      rsquared[i, j] <- s$r.squared[1]
    }
  }
  prop <- (pr$sdev[1:npc]^2/sum(pr$sdev^2)) * 100
  assoc <- linp
  colnames(assoc) <- paste(rep("PC", npc), 1:npc, " (", round(prop[1:npc],2), "%) ", sep="")
  return(assoc)
}

assoc.var.int <- function(X, pd, verbose = TRUE, ...) {
  # run the principal-component analysis
  var.names <- colnames(pd) 
  # remove variables with one level
  test <- sapply(var.names, function(v) { 
    if(is.factor(pd[,v]))
      if(nlevels(pd[,v]) == 1) TRUE
    else FALSE
    else FALSE}) 
  ids <- which(test == TRUE)
  if(length(ids) > 0) var.names <- var.names[-ids]
  # compile assoc matrix
  assoc <- matrix(NA, nrow=length(var.names), ncol=ncol(X))
  rownames(assoc) <- var.names
  colnames(assoc) <- colnames(colnames(X))
  print("assoc:")
  print(assoc)
  print("str(assoc):")
  print(str(assoc))
  for(i in 1:ncol(assoc)) {
    cat("assoc col : ", i, "\n")
    assoc[,i] <- sapply(var.names, function(a) {
      cat("pd var name : ", a, "\n")
      if(is.factor(pd[,a])) {
        cat("Is a factor\n")
        if(length(levels(pd[,a])) > 0)
          anova(lm(X[,i]~as.factor(pd[,a])))[1,5]
      }
      else if(is.numeric(pd[,a]))
        cat("Is numeric\n")
      anova(lm(X[,i]~pd[,a]))[1,5]
    })
  }
  return(assoc)
}

remove.batch.effects <- function(comb.data, pd, num.pc, vars, inter = "", method = "Combat", plot = TRUE, verbose = TRUE, ...) {
  batches <- vars$batches
  covs <- vars$covariates
  if(plot){
    cat(" - Running Principal Component Analysis..\n")
    assoc <- pc.anaylsis.2(g = comb.data, o = pd, npc = num.pc)
    print(assoc)
    limma:::plotMDS(comb.data, top=500, labels=pd[,vars$var.int], col=as.numeric(pd[,vars$var.int]), gene.selection="common", main = "Before removing any batch.")
  }
  for(batch in batches){
    id.var <- grep(batch, batches)[1]
    print(id.var)
    pre.batches <- batches
    batches <- batches[-id.var]
    string.formula <- paste("~", inter, "pd$", vars$var.int, sep="")
    if(length(batches) > 0)
      string.formula <- paste(string.formula, "+pd$", paste(batches, collapse="+pd$"), sep="")
    if(length(covs) > 0)
      string.formula <- paste(string.formula, "+pd$", paste(covs, collapse="+pd$"), sep="")
    form <- formula(string.formula)
    print(form)
    prev.comb.data <- comb.data
    if(method == "Combat") 
      res <- try(comb.data <- sva::ComBat(dat=comb.data, batch=pd[,batch], 
                                          mod=model.matrix(form), 
                                          par.prior=T,  prior.plots=F))
    if(method == "Limma")
      res <- try(comb.data <- limma::removeBatchEffect(x = comb.data,
                                                       batch = pd[,batch],
                                                       design = model.matrix(form)))
    
    if(inherits(res, "try-error")) {
      ##error handling: restore the previous state
      #comb.data <- prev.comb.data 
      #batches <- pre.batches
      #print(c(vars$var.int,batches))
      #next
      
      #error handling: return the custom error statement, with variable names
      error.str <- paste0("Variables are confounded!\n\nPlease remove one or more variable, so the design is not confunded!\n\nError encountered while correcting for BATCH:", batch, "\n\nPlease check ComBat resource for more information!")
      return(error.str)
    }
    if(plot){
      cat(" - Running Principal Component Analysis..\n")
      assoc <- pc.anaylsis.2(g=comb.data, o=pd, npc=num.pc)
      print(assoc)
      limma:::plotMDS(comb.data, top=500, labels=pd[,vars$var.int], col=as.numeric(pd[,vars$var.int]),gene.selection="common", main = paste("After removing", batch))
    }
  }
  if(plot)  {
    boxplot(comb.data, las=2, cex=0.7, main="After removing batch effects.")
    #plot.heatmap(comb.data, pd[,vars$var.int], 1000)
  }
  return(comb.data)
}

remove.batch.effects.2 <- function(comb.data, pd, num.pc, vars, inter = "", method = "Combat", plot = TRUE, verbose = TRUE, ...) {
  batches <- vars$batches
  covs <- vars$covariates
  cat(" - Running Principal Component Analysis..\n")
  assoc <- pc.anaylsis(data = comb.data, pd = pd, npc = num.pc)
  print(assoc)
  if(plot) limma:::plotMDS(comb.data, top=1000, labels=pd[,vars$var.int], col=as.numeric(pd[,vars$var.int]), 
                           gene.selection="common", main = "Before removing any batch.")
  #limma:::plotMDS(comb.data, top=500, labels=pd[,"TimeExposition"], col=as.numeric(pd[,"TimeExposition"]), 
  #                gene.selection="common", main = "Before removing any batch.")
  flag <- continue.read.inputs()
  while(flag) {
    batch  <- read.input.word(message = " - Please enter a variable <name> to use as batch..\n", 
                              more.info = paste("batches:", paste(batches, collapse=",")))
    id.var <- which(batches == batch)
    print(id.var)
    if(length(id.var) > 0) {
      prev.comb.data <- comb.data
      pre.batches <- batches
      print(c(covs, batches)[(length(covs) + id.var)])
      comb.data <- swamp:::combat(comb.data, pd[,c(covs, batches)], length(covs) + id.var, par.prior=T,  prior.plots=F)
      batches <- batches[-id.var]
      cat(" - Running Principal Component Analysis..\n")
      assoc <- pc.anaylsis(data = comb.data, pd = pd, npc = num.pc)
      print(assoc)
      if(plot) limma:::plotMDS(comb.data, top=1000, labels=pd[,vars$var.int], col=as.numeric(pd[,vars$var.int]), 
                               gene.selection="common", main = paste("After removing", batch))
      # roll back if the user wants to
      rb <- continue.read.inputs(message = "Do you want to roll back? (y/n)")
      if(rb == TRUE) {
        comb.data <- prev.comb.data 
        batches <- pre.batches
        cat(" - Previous Principal Component Analysis..\n")
        assoc <- pc.anaylsis(data = comb.data, pd = pd, npc = num.pc)
        print(assoc)
      }
    }
    else {
      cat(" - Invalid variable name.\n")
      next
    }
    flag <- continue.read.inputs()
  }
  if(plot)  boxplot(comb.data, las=2, cex=0.7, main="After removing batch effects.")
  return(comb.data)
}

get.sva.batch.effects <- function(comb.data, pd, vars, npc = 10, verbose = T, cmd.ui = T) {
  print("Searching for unknown batch effects:")
  covs <- vars$covariates
  string.formula <- paste("~pd$", vars$var.int, sep="")
  if(length(covs) > 0)  string.formula <- paste(string.formula, "+pd$", paste(covs, collapse="+pd$"), sep="")
  print(string.formula)
  form <- formula(string.formula)
  print(form)
  X <- sva::sva(dat=comb.data, mod=model.matrix(form), method = "two-step")$sv
  cat("class(X) - ", class(X), "\n")
  cat("dim(X) - ", dim(X), "\n")
  X.c <- X
  cat("class(X.c) - ", class(X.c), "\n")
  cat("dim(X.c) - ", dim(X.c), "\n")
  X <- infotheo::discretize(as.matrix(X), disc="equalfreq", nbins=NROW(X)^(1/3))
  X <- as.data.frame(X)
  if(verbose) print(X)
  colnames(X) <- paste("sva",c(1:ncol(X)),sep=".")
  #colnames(X) <- paste("svaD",c(1:ncol(X)),sep=".")
  rownames(X) <- colnames(comb.data)
  #colnames(X.c) <- paste("svaC",c(1:ncol(X)),sep=".")
  colnames(X.c) <- colnames(X)
  rownames(X.c) <- rownames(X)
  assoc.pca <- pc.anaylsis.2(g = comb.data, o = X, npc = npc)
  if(verbose) print(assoc.pca)
  assoc.pd <- assoc.var.int(X, pd)
  colnames(assoc.pd) <- colnames(X)
  if(verbose) print(assoc.pd)
  # indicating the batch variables found by uysing SVA
  if(cmd.ui){
    batches <- read.input.word(message = " - indicate batches <n to exit>: ")
    if(batches != "n") {
      ids <- as.numeric(unlist(strsplit(batches, ",")))
      X[,ids]
    }else NULL
  }else list(pd=assoc.pd, sv=X, svc=X.c)
}

## ANNOT-PROCESSING functions
connect.to.biomart <- function(verbose = TRUE, ...) {
  marts <- listMarts(host="www.ensembl.org")
  print(as.data.frame(marts[,1]))
  n1 <- read.input.number("Please enter a number to indicate the BioMart database to use.")
  if(verbose) cat(" - selected BioMart database:", as.character(marts[n1,1]), "\n", sep="")
  mart <- useMart(as.character(marts[n1,1]))
  datasets <- listDatasets(mart)
  organism <- read.input.word(1, "mmusculus or hsapiens")
  print(class(organism))
  ids <- grep(organism, datasets[,1])
  if(length(ids) > 0) cat(" - selected organism:", as.character(datasets[ids,1]), "\n", sep="")
  else stop("dataset not found")
  mart <- useMart(biomart=as.character(marts[n1,1]), dataset=as.character(datasets[ids,1]))
  mart
}

connect.to.biomart.2 <- function(verbose = TRUE, ...) {
  organisms <- c("hsapiens_gene_ensembl","mmusculus_gene_ensembl") 
  print(as.data.frame(organisms, stringsAsFactors = F))
  cat("Please select the BioMart database to use...\n")
  n1 <- scan(n=1)     
  if(is.numeric(n1)) 
    list(mart = useMart("ENSEMBL_MART_ENSEMBL", dataset=organisms[n1], host="www.ensembl.org"), organisms[n1])
  else
    NULL
}

query.biomart <- function(mart, organism, data, verbose = TRUE, ...) {
  probes <- as.character(unique(rownames(data)))
  if(interactive()){
    agilent.annot <- c("efg_agilent_wholegenome_4x44k_v1",  "efg_agilent_wholegenome_4x44k_v2",
                       "efg_agilent_sureprint_g3_ge_8x60k", "efg_agilent_sureprint_g3_ge_8x60k_v2")
    print(as.data.frame(agilent.annot, stringsAsFactors = F))
    cat("Please enter a number to indicate the filter to use as input to the query...\n")
    n1 <- scan(n=1)
    attributes <- listAttributes(mart)[c(1:100),1]
    print(as.data.frame(attributes, stringsAsFactors = F))
    cat("Please enter a number to indicate the annotation to use...\n")
    n2 <- scan(n=1)
    resq <- NULL       
    if(is.numeric(n1) & is.numeric(n2)) {
      # querying BioMart
      print(c(agilent.annot[n1], attributes[n2]))
      print(agilent.annot[n1])
      print(head(probes))
      resq <- getBM(attributes = c(agilent.annot[n1], attributes[n2]), 
                    filters = agilent.annot[n1], 
                    values = probes, 
                    mart = mart)
    }
    else stop("Invalid inputs.")
    # keep only complete cases
    if(verbose) cat(" - number of probes referring to 'NA' values:", length(which(complete.cases(resq) == FALSE)),"\n", sep="")
    if(length(which(complete.cases(resq) == FALSE)) > 0) resq <- resq[complete.cases(resq),]
    # remove probes referring to different gene ids
    u.probes <- unique(resq[,1])
    tests <- unlist(lapply(u.probes, function(x) { length(resq[which(x == resq[,1]), 2]) > 1}))
    if(verbose) cat(" - number of shared probes:", length(which(tests == TRUE)),"\n", sep="")
    if(length(which(tests == TRUE)) > 0) {
      u.probes <- u.probes[which(tests == FALSE)]
      ids.to.keep <- unlist(lapply(u.probes, function(p) which(resq[,1] == p)))
      resq <- resq[ids.to.keep,]
    }
    if(verbose) cat(" - number of unique probes:", length(unique(resq[,1])), "\n", sep="")
    # check species
    testh <- grep("hsap", organism)
    if(length(testh) > 0) spec <- "Homo sapiens"
    else spec <- "##"
    print(head(resq[,2]))
    print(attributes[n2])
    annotations <- get.genes(geneids = resq[,2], mart = mart, species=spec, type = attributes[n2])
    return(list(map = resq, annot = annotations))
  }
  NULL
}

getGene <- function (id, type, mart) {
  #biomaRt:::martCheck(mart, "ensembl")
  biomaRt:::checkWrapperArgs(id, type, mart)
  symbolAttrib = switch(strsplit(biomaRt:::martDataset(mart), "_", fixed = TRUE, 
                                 useBytes = TRUE)[[1]][1], hsapiens = "hgnc_symbol", mmusculus = "mgi_symbol", "external_gene_id")
  typeAttrib = switch(type, affy_hg_u133a_2 = "affy_hg_u133a_v2", type)
  attrib = c(typeAttrib, symbolAttrib, "description", "chromosome_name", 
             "band", "strand", "start_position", "end_position", "ensembl_gene_id")
  table = biomaRt:::getBM(attributes = attrib, filters = type, values = id, mart = mart)
  return(table)
}

get.genes <- function(geneids, mart, species="Homo sapiens", type) {
  require(biomaRt)
  if(species=="Homo sapiens")
    gannot <- getGene(id = geneids, type = type, mart = mart)[,c(type, "hgnc_symbol","description","ensembl_gene_id")]
  else {
    gannot <- getGene(id = geneids, type = type, mart = mart)  #[,c(type, "symbol","description","ensembl_gene_id")]
    print(colnames(gannot))
  }
  id.dups <- which(duplicated(gannot[,type]))
  if(length(id.dups) > 0) gannot <- gannot[-id.dups,]
  print(head(gannot)) 
  gannot[,"description"] <- gsub(" \\[.*\\]", "", gannot[,"description"])
  rownames(gannot) <- gannot[,type]
  gannot
}

## GENE EXPRESSION ANALYSIS
aggreg.probes <- function(data, map, var.mds = NULL, plot = TRUE, verbose = TRUE, ...) {
  genes <- unique(map[,2])
  gene.to.probes <- lapply(genes, function(g) map[which(map[,2] == g),1])
  names(gene.to.probes) <- genes
  probes.to.dupl <- lapply(unlist(gene.to.probes), function(probe) which(rownames(data) == probe))
  names(probes.to.dupl) <- unlist(gene.to.probes)
  mat <- apply(data, 2, function(s) {
    sapply(gene.to.probes, 
           function(probes) {
             #print(unlist(probes.to.dupl[probes]))
             median(s[unlist(probes.to.dupl[probes])])
           })
  })
  if(verbose) cat("Number of annotated genes: ", nrow(mat), "\n", sep="")
  if(plot & !is.null(var.mds)) 
    limma:::plotMDS(mat, top=1000, labels=var.mds, col=as.numeric(var.mds), 
                    gene.selection="common", main = "After mapping probes to genes")
  mat
}

aggreg.probes.2 <- function(data, map, var.mds = NULL, plot = TRUE, verbose = TRUE, ...) {
  #Update map column 2 to TargetID
  colnames(map)[2] <- "TargetID"
  #Aggregate dup probes before aggregating by annotation IDs
  tmpRowNames <- rownames(data)
  dups <- which(duplicated(tmpRowNames))
  if(length(dups)>0){
    dupRowNames <- unique(tmpRowNames[dups])
    unique.data <- data[-which(tmpRowNames %in% dupRowNames),]
    colnames(unique.data) <- make.names(colnames(unique.data))
    dup.data <- data[which(tmpRowNames %in% dupRowNames),]
    dup.data <- data.frame(ProbeID=rownames(dup.data), dup.data, row.names=NULL)
    dup.datag <- aggregate(.~ProbeID, data = dup.data, median)
    rownames(dup.datag) <- dup.datag[,1]
    dup.datag <- dup.datag[,-1]
    data <- rbind(unique.data, dup.datag)
  }
  rownames(map) <- map[,1]
  map <- map[rownames(data),]
  map <- cbind(map,data)
  map <- map[,-1]
  blankIdx <- which(map$TargetID=="")
  if(length(blankIdx)>0){
    map <- map[-blankIdx,]
  }
  if(verbose)  print(str(map))
  datag = aggregate(.~TargetID, data = map, median)
  rownames(datag) <- datag[,1]
  datag <- datag[,-1]
  if(plot & !is.null(var.mds))  
    limma:::plotMDS(datag, top=1000, labels=var.mds, col=as.numeric(var.mds), gene.selection="common", main = "After mapping probes to genes")
  datag
}

diff.gene.expr <- function(data, des, contrasts, pvalue, fcvalue, p.adjust.method, annot = NULL, max.ngenes=Inf, save.file = NULL, plot = TRUE, verbose = TRUE, ...) {
  # make the matrix of contrasts
  colnames(des) <- make.names(colnames(des))
  cont <- limma::makeContrasts(contrasts=contrasts, levels=des)
  if(verbose) {
    print(des)
    print(cont)
  }
  fit <- limma::lmFit(data, des)
  ##Added the print for contrast matrix
  #print(contrasts.fit(fit, cont))
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont))
  #max.ngenes <- read.input.number(" - indicate the max number of genes:")
  if(plot & length(contrasts) < 6) {
    limma::vennDiagram(limma::decideTests(fit2, p.value=pvalue, lfc=log2(fcvalue), adjust.method=p.adjust.method),
                       cex=c(.8,.8,0.7),
                       main = paste("Venn Diagram <p.value=", pvalue, ",FC=", fcvalue, ">",sep=""))
    #vennDiagram(decideTests(fit2, p.value=pvalue, lfc=log2(fcvalue), adjust.method=methods[n1]), include=c("up","down"), counts.col=c("red","green"),
    #            main = paste("Venn Diagram <p.value=", pvalue, ",FC=", fcvalue, ">",sep=""))
    ##for(i in 1:length(vars$contrasts))  volcanoplot(fit2, coef=i, highlight=50)
  }
  # to output a table of the differentially expressed genes
  list.top.tables <- list()
  if(verbose)  cat(" - generating top tables.. \n", sep="")
  for(i in 1:length(contrasts)) {
    list.top.tables[[i]] <- limma::topTable(fit2, coef=i, p.value=pvalue, lfc=log2(fcvalue), 
                                            adjust.method=p.adjust.method,
                                            number=max.ngenes, sort.by="logFC")
    if(verbose) {
      cat("     > ", contrasts[i], " - number of sig. genes:",  nrow(list.top.tables[[i]]), "\n", sep="")
    }
  }
  names(list.top.tables) <- contrasts
  
  # how to order the gene expression data
  for(i in 1:length(list.top.tables)) {
    ids <- NULL
    if(dim(list.top.tables[[i]])[1] == 0) next
    if(!is.null(annot)) {
      list.top.tables[[i]] <- base:::merge(list.top.tables[[i]], annot, by = "row.names")
      rownames(list.top.tables[[i]]) <- list.top.tables[[i]]$Row.names
      list.top.tables[[i]] <- list.top.tables[[i]][,-1]
      #list.top.tables[[i]]$score <- abs(list.top.tables[[i]]$logFC) * -log(list.top.tables[[i]]$P.Value)
    }else{
      ids <- rownames(list.top.tables[[i]])
    }
    list.top.tables[[i]]$score <- list.top.tables[[i]]$logFC * -log10(list.top.tables[[i]]$P.Value)
    if(!is.null(ids)) list.top.tables[[i]]$ID <- rownames(list.top.tables[[i]])
    colnames(list.top.tables[[i]])[which(colnames(list.top.tables[[i]]) == "t")] <- "t-statistic" 
    colnames(list.top.tables[[i]])[which(colnames(list.top.tables[[i]]) == "B")] <- "B-statistic"
  }
  list.top.tables
}

volcano.plot <- function(data, cutoff.pvalue = .05, cutoff.lfc = 1) {
  require(calibrate)
  plot(data$logFC, -log10(data$P.Value), pch=20, main="Volcano plot", xlim=c(-2.5,2))
  pdata <- subset(data, P.Value < cutoff.pvalue) 
  points(pdata$logFC, -log10(pdata$P.Value), pch=20, col="gray")
  fdata <- subset(data, abs(logFC) > cutoff.lfc)
  points(fdata$logFC, -log10(fdata$P.Value), pch=20, col="orange")
  gdata <- subset(data, P.Value < cutoff.pvalue & abs(logFC) > cutoff.lfc)
  points(gdata$logFC, -log10(gdata$P.Value), pch=20, col="blue")
  textxy(gdata$logFC, -log10(gdata$P.Value), labs=gdata$hgnc_symbol, cex=.6)
}

save.xls.file <- function(lists.diff.genes = NULL, file.name = "diff_genes.xlsx") {
  #require(xlsx)
  require(XLConnect)
  options(java.parameters = "-Xmx4g")
  if(!is.null(lists.diff.genes)) {
    gc()
    # creating work book
    # wb <- createWorkbook()
    xlcFreeMemory()
    wb <- loadWorkbook(file.name, create=TRUE)
    for(i in 1:length(lists.diff.genes)) {
      # making each as a sheet
      #sheet <- createSheet(wb, sheetName=names(lists.diff.genes)[i])
      sheet <- names(lists.diff.genes)[i]
      createSheet(wb, sheet)
      #addDataFrame(lists.diff.genes[[i]], sheet)
      writeWorksheet(wb, lists.diff.genes[[i]], sheet, startRow=1, startCol=1, header=TRUE, rownames=rownames(lists.diff.genes[[i]]))
    }
    # saving the workbook
    #saveWorkbook(wb, file.name)
    saveWorkbook(wb)
  }
  NULL
}

david.annot <- function(list.top.tables) {
  require(BACA)
  require(ggplot2)
  split <- read.input.number("Enter [1] to split each list in two sublists: down- and up-regulated genes.")
  # building the gene list
  if(split == 1) {
    gene.lists <- unlist(lapply(list.top.tables, function(d) { 
      down.genes <- rownames(d)[which(d$logFC < 0)]
      if(length(down.genes)==0) down.genes <- "1"
      up.genes <- rownames(d)[which(d$logFC >= 0)]  
      if(length(up.genes)==0) up.genes <- "1"
      list(down.genes, up.genes)}), recursive = FALSE)
    # naming the gene lists 
    names(gene.lists) <- unlist(lapply(names(list.top.tables), function(x) paste(x,c(1,2),sep="_")))
    print(length(gene.lists))
  }
  else {
    gene.lists <- lapply(list.top.tables, rownames)
    names(gene.lists) <- names(list.top.tables)
  }
  return(gene.lists)
  # querying DAVID for KEGG.PATHWAYS
  kegg.pathways <- DAVIDsearch(gene.lists, david.user = "vittorio.fortino@ttl.fi", 
                               idType = "ENTREZ_GENE_ID", annotation = "KEGG_PATHWAY")
  
  bp.terms <- DAVIDsearch(gene.lists, david.user = "vittorio.fortino@ttl.fi", 
                          idType = "ENTREZ_GENE_ID", annotation = "GOTERM_BP_ALL")
  
  mf.terms <- DAVIDsearch(gene.lists, david.user = "vittorio.fortino@ttl.fi", 
                          idType = "ENTREZ_GENE_ID", annotation = "GOTERM_MF_ALL")
  
  cc.terms <- DAVIDsearch(gene.lists, david.user = "vittorio.fortino@ttl.fi", 
                          idType = "ENTREZ_GENE_ID", annotation = "GOTERM_CC_ALL")
  
  list(show.bbplot(kegg.pathways, "BBplot - KEGG", names(list.top.tables)),
       show.bbplot(bp.terms, "BBplot - GO terms (BP)", names(list.top.tables)),
       show.bbplot(mf.terms, "BBplot - GO terms (MF)", names(list.top.tables)),
       show.bbplot(cc.terms, "BBplot - GO terms (CC)", names(list.top.tables)))
  
} 

show.bbplot <- function(res.david.query, title.query, col.names, verbose = TRUE) {
  cont.it <- TRUE
  if(verbose)  cat("Building the bbplot..", title.query, "\n", sep="")
  while(cont.it) {
    max.pval <- read.input.number("Indicate the p.value:")
    min.ngenes <- read.input.number("Indicate the min number of genes:")
    cat("Do you want to specify the max number of genes?", "\n")
    if(continue.read.inputs()) {
      max.ngenes <- read.input.number("Indicate the max number of genes:")
      bbplot <- BBplot(res.david.query, max.pval = max.pval, 
                       min.ngenes = min.ngenes, 
                       max.ngenes = max.ngenes,
                       adj.method = "",
                       name.com = col.names, 
                       title = paste(title.query, "<p.value=", max.pval,
                                     ";min.num.genes=", min.ngenes,
                                     ";max.num.genes=", max.ngenes,
                                     ">", sep=""), print.term = "full") 
      print(bbplot)
    }
    else {
      # plotting
      bbplot <- BBplot(res.david.query, max.pval = max.pval, min.ngenes = min.ngenes, name.com = col.names, 
                       title = paste(title.query, "<p.value=", max.pval,";min.num.genes=", min.ngenes,">", sep=""), print.term = "full") 
      print(bbplot)
    }
    cont.it <- continue.read.inputs()
  }
  bbplot
}

## USER INPUTS FUNCTIONS
read.input.word <- function(n=1, message = "Indicate a number", more.info=NULL) {
  cat(message, "\n", sep="")
  if(!is.null(more.info)) cat(" *** ", more.info, "\n")
  keywords <- scan(what=character(), nlines=n)
  paste(keywords, collapse=",")
}

read.input.number <- function(message = "Indicate a number", more.info=NULL) {
  cat(message, "\n", sep="")
  if(!is.null(more.info)) cat(more.info, "\n", sep="")
  n <- NULL
  while(!is.numeric(n)) {
    n <- scan(n=1)
    if(!is.numeric(n))
      cat("Invalid input number.", "\n")
  }
  n
}

continue.read.inputs <- function(message="Press [y] to continue, [n].") {
  next.it <- ""
  while(next.it != "y" &  next.it != "n") {
    cat(message)
    next.it <- read.input.word(message = "")
  }
  if(next.it == "y") TRUE
  else FALSE
} 

## UTILITY functions
droplevels.annot <- function(a, remove = F) {
  id.col.to.remove <- NULL
  for(j in 1:ncol(a)) {
    if(is.character(a[,j]))
      a[,j] <- as.factor(a[,j])
    if(is.factor(a[,j])) {
      a[,j] <- droplevels(a[,j])
    }
    if(length(levels(a[,j])) > 1)
      print(colnames(a)[j])
    else
      id.col.to.remove <- c(id.col.to.remove, j)
  }
  print(j)
  if(remove) 
    a[,-j]
  else
    a
}

build.model.matrix <- function(pd, intercept = -1, var.int, covariates = NULL, verbose = TRUE) {
  des <- NULL
  if(!is.null(covariates)) {
    f <- formula(paste("~", intercept, "+pd$", var.int, "+pd$", paste(covariates, collapse="+pd$"), sep=""))
    if(verbose) print(f)
    des <- model.matrix(f)
    colnames(des) <- gsub(paste("pd\\$", var.int, sep=""), "", colnames(des))
    for(i in 1:length(covariates))
      colnames(des) <- gsub( paste("pd\\$", covariates[i], sep=""), "", colnames(des))
  }
  else {
    f <- formula(paste("~", intercept, "+pd$", var.int, sep=""))
    if(verbose) print(f)
    des <- model.matrix(f)
    colnames(des) <- levels(as.factor(pd[,var.int]))
  }
  des
}

color.map <- function(labels) {
  levs <- levels(labels)
  colors <- rep(NA, length(labels))
  for(i in 1:length(levs)) {
    colors[which(labels == levs[i])] <- COLORS[i]
  }
  colors
}

hclust2 <- function(x, method="ward.D2", ...) hclust(x, method=method, ...)
dist2 <- function(x, ...) as.dist(1-cor(t(x), method="pearson"))
plot.heatmap <- function(data, labels, top = 500) {
  sidebarcolors <- color.map(labels)
  hvars <- apply(data, 1, var,na.rm=TRUE) 
  hvars <- sort(hvars, decreasing=TRUE) 
  hvars <- hvars[1:top] 
  print(dim(data[names(hvars),]))
  heatmap.2(data[names(hvars),], Rowv=TRUE, 
            scale="column", trace="none", 
            distfun=dist2, hclustfun = hclust2,
            col=redgreen, xlab=NULL, ylab=NULL, 
            labRow = "", dendrogram = "column", cexCol = .7,
            ColSideColors=sidebarcolors)
  par(lend = 1)           # square line ends for the color legend
  legend("bottomleft",      # location of the legend on the heatmap plot
         legend = levels(labels), # category labels
         col = COLORS[1:length(levels(labels))],  # color key
         bty = "n",
         bg = "white",
         border = "white",
         cex = .7,
         lty= 1,            # line style
         lwd = 10           # line width
  )     
}

plot.heatmpa.genes <- function(expr.dat, col, sf = 3, ze = 0.4, pi = 0.4, o = 1, ccolors) {
  require(Biobase)
  require(DFP)
  dfp <- featSelDFP(expr.dat, skipFactor = sf, zeta = ze, piVal = pi, overlapping = o)[[3]]
  dl <- as.matrix(dfp)
  dl <- apply(dl, MARGIN = 2, FUN = function(x) {
    id.na <- which(is.na(x))
    if(length(id.na) > 0) x[id.na] <- "Not assigned"
    x
  })
  feats <- rownames(dfp)
  dat <- expand.grid(var1=1:nrow(dfp), var2=1:ncol(dfp))
  lev <- colnames(attr(dfp,"ifs"))
  dat$var2 <- sapply(dat$var2, function(i) lev[i])
  dat$value <- round(melt(attr(dfp,"ifs"))$value, 2)
  dat$labels <- melt(dl)$value
  print(dat)
  ggplot(dat, aes(x=var1,y=var2))+ 
    geom_point(aes(size = value, colour = labels)) +
    scale_colour_manual(values=ccolors) +
    scale_size(range = c(5, 10), breaks=c(0.25,.5,.75),labels=c("25%","50%","75%")) +
    #scale_size(range = c(5, 10)) + 
    scale_x_continuous(breaks = 1:nrow(dfp), labels = feats) +
    scale_y_discrete(limits = sort(levels(expr.dat$class), decreasing = T)) + 
    labs(size = "Coverage", colour="Gene status") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          panel.background = element_blank(),
          legend.text = element_text(size = 10),
          legend.position="right") +
    xlab("Genes") +
    ylab("Class Labels") 
}

split.top.tables <- function(list.tables, fc = 1.5) { 
  gene.lists <- unlist(lapply(list.tables, function(d) { 
    rownames(d) <- sub("_at",'', rownames(d))
    down.genes <- rownames(d)[which(d$logFC < -log2(fc))]
    if(length(down.genes)==0) down.genes <- "1"
    up.genes <- rownames(d)[which(d$logFC >= log2(fc))]  
    if(length(up.genes)==0) up.genes <- "1"
    list(down.genes, up.genes)}), recursive = FALSE)
  # naming the gene lists 
  names(gene.lists) <- unlist(lapply(names(list.tables), function(x) paste(x,c("down","up"),sep="_")))
  print(length(gene.lists))
  print(lapply(gene.lists, length))
  gene.lists
}

# Goodman and Kruskal's tau measure
GKtau <- function(x,y){
  #  First, compute the IxJ contingency table between x and y
  Nij = table(x,y,useNA="ifany")
  #  Next, convert this table into a joint probability estimate
  PIij = Nij/sum(Nij)
  #  Compute the marginal probability estimates
  PIiPlus = apply(PIij,MARGIN=1,sum)
  PIPlusj = apply(PIij,MARGIN=2,sum)
  #  Compute the marginal variation of y
  Vy = 1 - sum(PIPlusj^2)
  #  Compute the expected conditional variation of y given x
  InnerSum = apply(PIij^2,MARGIN=1,sum)
  VyBarx = 1 - sum(InnerSum/PIiPlus)
  #  Compute and return Goodman and Kruskal's tau measure
  tau = (Vy - VyBarx)/Vy
  tau
}

# DFP-based analysis
featSelDFP <- function(exprData, skipFactor = 3, zeta = 0.5, piVal = 0.9, overlapping = 2) {
  require(DFP)
  mfs <- calculateMembershipFunctions(exprData, skipFactor)
  print("calculateMembershipFunctions <DONE>")
  dvs <- discretizeExpressionValues(exprData, mfs, zeta, overlapping)
  print("discretizeExpressionValues <DONE>")
  fps <- calculateFuzzyPatterns(exprData, dvs, piVal, overlapping)
  print("calculateFuzzyPatterns <DONE>")
  list.fps <- lapply(names(table(exprData$class)), function(c) { 
    fp <- showFuzzyPatterns(fps, c) 
    rownames(exprData)[which((rownames(exprData) %in% names(fp)) == TRUE)]
  })
  #xlistFPs <- lapply(names(table(exprData$class)), FUN= function(c) { fp = showFuzzyPatterns(fps, c); which( (rownames(exprData) %in% names(fp)) == TRUE)})
  # print(xlistFPs)
  dfps <- calculateDiscriminantFuzzyPattern(exprData, fps)
  #plotDiscriminantFuzzyPattern(dfps, overlapping)
  return(list(fps,list.fps,dfps,dvs))
}

merge.combat <- function (esets, covariates = NULL, verbose = T)  {
  raw_merged = merge.array(esets)
  batchInfo = NULL
  for (i in 1:length(esets)) {
    batchInfo = c(batchInfo, rep(i, ncol(esets[[i]])))
  }
  #mod = cbind(as.factor(imp_ann$Nano2), as.factor(imp_ann$TimeExposition),as.factor(imp_ann$cell_type2))
  #print(colnames(pData(raw_merged)))
  c.names <- c("Array name", "Sample name", "Batch")
  saminfo = cbind(rownames(pData(raw_merged)), rownames(pData(raw_merged)), batchInfo)
  if(!is.null(covariates)) {
    for(i in 1:length(covariates)) {
      saminfo <- cbind(saminfo, as.factor(pData(raw_merged)[,covariates[i]]))
      c.names <- c(c.names, covariates[i])
    }
  }
  #colnames(saminfo) = c("Array name", "Sample name", "Batch")
  colnames(saminfo) = c.names
  dat = exprs(raw_merged)
  design <- design.mat(saminfo)
  if(verbose) print(design)
  batches <- list.batch(saminfo)
  n.batch <- length(batches)
  n.batches <- sapply(batches, length)
  n.array <- sum(n.batches)
  B.hat <- solve(t(design) %*% design) %*% t(design) %*% t(as.matrix(dat))
  grand.mean <- t(n.batches/n.array) %*% B.hat[1:n.batch, ]
  var.pooled <- ((dat - t(design %*% B.hat))^2) %*% rep(1/n.array, 
                                                        n.array)
  stand.mean <- t(grand.mean) %*% t(rep(1, n.array))
  if (!is.null(design)) {
    tmp <- design
    tmp[, c(1:n.batch)] <- 0
    stand.mean <- stand.mean + t(tmp %*% B.hat)
  }
  s.data <- (dat - stand.mean)/(sqrt(var.pooled) %*% t(rep(1, n.array)))
  batch.design <- design[, 1:n.batch]
  gamma.hat <- solve(t(batch.design) %*% batch.design) %*% 
    t(batch.design) %*% t(as.matrix(s.data))
  delta.hat <- NULL
  for (i in batches) {
    delta.hat <- rbind(delta.hat, apply(s.data[, i], 1, var, 
                                        na.rm = T))
  }
  gamma.bar <- apply(gamma.hat, 1, mean)
  t2 <- apply(gamma.hat, 1, var)
  a.prior <- apply(delta.hat, 1, aprior)
  b.prior <- apply(delta.hat, 1, bprior)
  gamma.star <- delta.star <- NULL
  for (i in 1:n.batch) {
    temp <- it.sol(s.data[, batches[[i]]], gamma.hat[i, ], 
                   delta.hat[i, ], gamma.bar[i], t2[i], a.prior[i], 
                   b.prior[i])
    gamma.star <- rbind(gamma.star, temp[1, ])
    delta.star <- rbind(delta.star, temp[2, ])
  }
  bayesdata <- s.data
  j <- 1
  for (i in batches) {
    bayesdata[, i] <- (bayesdata[, i] - t(batch.design[i, 
    ] %*% gamma.star))/(sqrt(delta.star[j, ]) %*% t(rep(1, 
                                                        n.batches[j])))
    j <- j + 1
  }
  bayesdata <- (bayesdata * (sqrt(var.pooled) %*% t(rep(1, 
                                                        n.array)))) + stand.mean
  eset = raw_merged
  exprs(eset) = bayesdata
  return(eset)
}

merge.array <- function (esets)  {
  eset1 = esets[[1]]
  annot1 = annotation(eset1)
  for (i in 2:length(esets)) {
    eset2 = esets[[i]]
    d1 = exprs(eset1)
    d2 = exprs(eset2)
    cg = sort(intersect(rownames(d1), rownames(d2)))
    if (length(cg) < (min(dim(d1)[1], dim(d2)[1])/100)) {
      msg(" ! WARNING ! Number of common genes < 1%")
    }
    fData = fData(eset1)[cg, ]
    p1 = pData(eset1)
    p2 = pData(eset2)
    cp = sort(intersect(colnames(p1), colnames(p2)))
    tp = sort(unique(union(colnames(p1), colnames(p2))))
    sp1 = setdiff(colnames(p1), cp)
    sp2 = setdiff(colnames(p2), cp)
    pheno = matrix(NA, ncol = length(tp), nrow = nrow(p1) + 
                     nrow(p2))
    rownames(pheno) = c(rownames(p1), rownames(p2))
    colnames(pheno) = tp
    if (length(cp) != 0) {
      pheno[1:nrow(p1), cp] = as.matrix(p1[, cp])
      pheno[(nrow(p1) + 1):(nrow(p1) + nrow(p2)), cp] = as.matrix(p2[, 
                                                                     cp])
    }
    if (length(sp1) != 0) {
      pheno[1:nrow(p1), sp1] = as.matrix(p1[, sp1])
    }
    if (length(sp2) != 0) {
      pheno[(nrow(p1) + 1):(nrow(p1) + nrow(p2)), sp2] = as.matrix(p2[, 
                                                                      sp2])
    }
    pData = as.data.frame(pheno)
    d1 = d1[cg, , drop = FALSE]
    d2 = d2[cg, , drop = FALSE]
    eset1 = new("ExpressionSet")
    exprs(eset1) = cbind(d1, d2)
    pData(eset1) = pData
    fData(eset1) = fData
    annot1 = c(annot1, annotation(eset2))
  }
  annotation(eset1) = unique(annot1)
  return(eset1)
}

design.mat <- function (saminfo, verbose = TRUE) {
  tmp <- which(colnames(saminfo) == "Batch")
  tmp1 <- as.factor(saminfo[, tmp])
  if(verbose) cat("  => Found", nlevels(tmp1), "batches")
  design <- build.design(tmp1, start = 1)
  ncov <- ncol(as.matrix(saminfo[, -c(1:2, tmp)]))
  if(verbose) cat("  => Found", ncov, "covariate(s)")
  if (ncov > 0) {
    for (j in 1:ncov) {
      tmp1 <- as.factor(as.matrix(saminfo[, -c(1:2, tmp)])[, 
                                                           j])
      design <- build.design(tmp1, des = design)
    }
  }
  design
}

build.design <- function (vec, des = NULL, start = 2) {
  tmp <- matrix(0, length(vec), nlevels(vec) - start + 1)
  for (i in 1:ncol(tmp)) {
    tmp[, i] <- vec == levels(vec)[i + start - 1]
  }
  cbind(des, tmp)
}

list.batch <- function (saminfo) {
  tmp1 <- as.factor(saminfo[, which(colnames(saminfo) == "Batch")])
  batches <- NULL
  for (i in 1:nlevels(tmp1)) {
    batches <- append(batches, list((1:length(tmp1))[tmp1 == levels(tmp1)[i]]))
  }
  batches
}

aprior <- function (gamma.hat){
  m = mean(gamma.hat)
  s2 = var(gamma.hat)
  (2 * s2 + m^2)/s2
}

bprior <- function (gamma.hat) {
  m = mean(gamma.hat)
  s2 = var(gamma.hat)
  (m * s2 + m^3)/s2
}

it.sol <- function (sdat, g.hat, d.hat, g.bar, t2, a, b, conv = 1e-04)  {
  n <- apply(!is.na(sdat), 1, sum)
  g.old <- g.hat
  d.old <- d.hat
  change <- 1
  count <- 0
  while (change > conv) {
    g.new <- postmean(g.hat, g.bar, n, d.old, t2)
    sum2 <- apply((sdat - g.new %*% t(rep(1, ncol(sdat))))^2, 
                  1, sum, na.rm = T)
    d.new <- postvar(sum2, n, a, b)
    change <- max(abs(g.new - g.old)/g.old, abs(d.new - d.old)/d.old)
    g.old <- g.new
    d.old <- d.new
    count <- count + 1
  }
  adjust <- rbind(g.new, d.new)
  rownames(adjust) <- c("g.star", "d.star")
  adjust
}

postmean <- function (g.hat, g.bar, n, d.star, t2)  {
  (t2 * n * g.hat + d.star * g.bar)/(t2 * n + d.star)
}

postvar <- function (sum2, n, a, b)  {
  (0.5 * sum2 + b)/(n/2 + a - 1)
}

## MDS plot
plotMDS <- function (eset, colLabel, symLabel, legend = TRUE, file = NULL,  ctr = FALSE, md = FALSE, ...) {
  if (!is.null(file)) {
    pdf(file, width = 12, height = 7)
  }
  mds = cmdscale(dist(t(exprs(eset))), eig = TRUE)
  colMap = makeColorMap(eset, colLabel)
  colVec = makeColorVec(eset, colLabel, colMap)
  tmp = par()$mar
  if (legend) {
    par(xpd = T, mar = par()$mar + c(0, 0, 0, 4))
  }
  range_x = range(mds$points[, 1])
  range_y = range(mds$points[, 2])
  plot(mds$points, col = colVec, pch = as.numeric(pData(eset)[, 
                                                              symLabel]), panel.first = {
                                                                U = par("usr")
                                                                rect(U[1], U[3], U[2], U[4], col = "azure2", border = "black", 
                                                                     lwd = 3)
                                                              }, lwd = 2, xlab = "", ylab = "", xlim = range_x, ylim = range_y, 
       ...)
  if(ctr) {
    rect(-.5, -.5, .5, .5, border = "black", lty = 3, lwd = 2)
    rect(-1.5, -1.5, 1.5, 1.5, border = "black", lty = 2, lwd = 2)
    rect(-2.5, -2.5, 2.5, 2.5, border = "black", lwd = 2)
  }
  if (legend) {
    x = range_x[2] + (range_x[2] - range_x[1]) * 0.05
    y = range_y[2] - (range_y[2] - range_y[1]) * 0.05
    syms = unique(pData(eset)[, symLabel])
    legend("topleft", legend = syms, pt.lwd = 2, cex = 0.7, pch = as.numeric(syms), 
           box.lwd = 3, bg = "azure2")
    legend("topright", inset=c(-0.3,0), legend = names(colMap), pt.lwd = 2, pch = 19, 
           col = unlist(colMap), box.lwd = 3, bg = "azure2")
    par(xpd = F, mar = tmp)
  }
  if (!is.null(file)) {
    dev.off()
  }
  if(md) return(mds)
}

makeColorMap <- function (eset, label){
  colMap = list()
  vec = unique(as.vector(pData(eset)[, label]))
  for (i in 1:length(vec)) {
    colMap[[vec[i]]] = COLORS[i]
  }
  return(colMap)
}

makeColorVec <- function (eset, label, colMap) {
  labels = as.vector(pData(eset)[, label])
  return(as.vector(unlist(sapply(labels, function(x) {
    colMap[x]
  }))))
}

COLORS <- c("green3","blue","cyan","magenta","yellow","gray","red","orange",
            "darkred", "green","darkblue","darkcyan","darkmagenta","darkorchid1",
            "darkgoldenrod3","aquamarine","darkslategray3","darkolivegreen3","lightcoral",
            "deeppink","gold","Olivedrab1","dimgrey","cornflowerblue","darkgreen",
            "burlywood3","steelblue4","orangered","purple1","khaki1","azure4",
            "blue1","blue2","blue3","blue4","coral1","coral2","coral3","coral4")

#Get Summary of Differential Expression Tables
get_deg_summary <- function(deg_list, names, lfc=0){
  cuts <- c(0.0001, 0.001, 0.01, 0.025, 0.05, 0.10, 1)
  cat("\nStatistical significance summary:\n")
  countsDF <- NULL
  deg_list <- as.list(deg_list)
  
  for(i in 1:length(deg_list)){
    objectName <- names(deg_list)[i]
    object <- deg_list[[objectName]]
    idx<-which(abs(object$logFC)>lfc); 
    if(length(idx)>0){
      object <- object[idx,]
      counts <- sapply(cuts, function(x) c("P.Value"=sum(object$P.Value < x), "adj.P.Val"=sum(object$adj.P.Val < x)))
    }else{
      counts <- data.frame(sapply(cuts, function(x){c(0,0)}))
      colnames(counts) <- cuts
      rownames(counts) <- c("P.Value", "adj.P.Val")
    }
    colnames(counts) <- paste("<", cuts, sep="")
    countName <- names[i]
    #rownames(counts) <- paste0(countName, ";", rownames(counts))
    ##print(counts)
    #countsDF <- rbind(countsDF, counts)
    print(class(counts))
    tmpDF <- data.frame(comparison=countName, score=rownames(counts), counts, stringsAsFactors=FALSE)
    colnames(tmpDF)[3:ncol(tmpDF)] <- paste("<", cuts, sep="")
    rownames(tmpDF) <- paste0(countName, ";", rownames(counts))
    countsDF <- rbind(countsDF, tmpDF)
    print("Summary Table Class:")
    print(class(countsDF))
  }
  return(countsDF)
}

#Get color palette for factor
get_color_palette <- function(iVec, asFactor=FALSE){
  set.seed(1) #Block randomness. Set seed
  print("Getting color palette..")
  if(asFactor){
    print("Input as factor...")
    if(is.factor(iVec)){
      print("Input is factor...")
      nm <- levels(iVec)
    }else if(is.vector(iVec)){
      print("Input is vector...")
      iVec <- factor(iVec)
      nm <- levels(iVec)
    }else{
      print("Neither factor nor vector, returning NULL!")
      return(NULL)
    }
    print("levels:")
    print(nm)
    #nmPalette <- setNames(randomcoloR::distinctColorPalette(length(nm)), seq(nm))
    nmPalette <- setNames(randomcoloR::randomColor(length(nm), luminosity="dark"), seq(nm))
    print("Palette base:")
    print(nmPalette)
    print("Levels as integer:")
    print(as.integer(iVec))
    colorVec <- nmPalette[as.integer(iVec)]
  }else{
    print("Input as vector...")
    if(is.factor(iVec)){
      print("Input is factor...")
      iVec <- as.character(iVec)
    }
    print("iVec:")
    print(iVec)
    #colorVec <- setNames(randomcoloR::distinctColorPalette(length(iVec)), iVec)
    colorVec <- setNames(randomcoloR::randomColor(length(iVec), luminosity="dark"), iVec)
  }
  print("Palette vector:")
  print(colorVec)
  return(colorVec)
}

#Function to get names of columns with datatype character
factorize_cols <- function(phTable, idx){
  for(i in idx){
    phTable = as.data.frame(phTable)
    phTable[, i] <- as.factor(as.vector(phTable[, i]))
  }
  return(phTable)
}

#holds the results of a QC comparison
setClass("QCStats",representation(scale.factors="numeric",target="numeric",percent.present="numeric",average.background="numeric",minimum.background="numeric",maximum.background="numeric",spikes="matrix",qc.probes="matrix",bioBCalls="character",arraytype="character"));

#accessor methods
setGeneric("sfs", function(object) standardGeneric("sfs"))
setMethod("sfs","QCStats",function(object) object@scale.factors)

setGeneric("target", function(object) standardGeneric("target"))
setMethod("target","QCStats",function(object) object@target)

setGeneric("percent.present", function(object) standardGeneric("percent.present"))
setMethod("percent.present","QCStats",function(object) object@percent.present)

setGeneric("avbg", function(object) standardGeneric("avbg"))
setMethod("avbg","QCStats",function(object) object@average.background)

setGeneric("minbg", function(object) standardGeneric("minbg"))
setMethod("minbg","QCStats",function(object) object@minimum.background)

setGeneric("maxbg", function(object) standardGeneric("maxbg"))
setMethod("maxbg","QCStats",function(object) object@maximum.background)

setGeneric("spikeInProbes", function(object) standardGeneric("spikeInProbes"))
setMethod("spikeInProbes","QCStats",function(object) object@spikes)

setGeneric("qcProbes", function(object) standardGeneric("qcProbes"))
setMethod("qcProbes","QCStats",function(object) object@qc.probes)

setGeneric("arrayType", function(object) standardGeneric("arrayType"))
setMethod("arrayType","QCStats",function(object) object@arraytype)


# alter the QC environment - effectively accessor functions into the environment

.qc.is.empty <- function() {
  get("empty",.qcEnv)
}

.qc.set.empty <- function(empty) {
  assign("empty",empty,envir=.qcEnv)
}


qc.set.array <- function(name) {
  .qc.set.empty(FALSE)
  assign("array",name,envir=.qcEnv)
}

qc.get.array <- function() {
  .qc.test()
  get("array",.qcEnv)
}

qc.get.tau <- function() {
  .qc.test()
  0.015;
}

qc.set.alpha1 <- function(value) {
  .qc.set.empty(FALSE)
  assign("alpha1",value,envir=.qcEnv)
}

qc.get.alpha1 <- function() {
  .qc.test()
  get("alpha1",.qcEnv)
}

qc.set.alpha2 <- function(value) {
  .qc.set.empty(FALSE)
  assign("alpha2",value,envir=.qcEnv)
}

qc.get.alpha2 <- function() {
  .qc.test()
  get("alpha2",.qcEnv)
}

qc.get.spikes <- function() {
  .qc.test()
  get("spikes",.qcEnv)
}

qc.get.spike <- function(name) {
  .qc.test()
  spikes <- get("spikes",.qcEnv)
  spikes[[name]]
}

qc.add.spike <- function(name,probeset) {
  .qc.set.empty(FALSE)
  spikes <- qc.get.spikes()
  names  <- names(spikes)
  if(name %in% names) {
    spikes[name] <- probeset
  }
  else {
    spikes <- c(spikes,probeset)
    names  <- c(names,name)
    names(spikes) <- names
  }
  assign("spikes",spikes,envir=.qcEnv)
}

qc.get.probes <- function() {
  .qc.test()
  get("probes",.qcEnv)
}

qc.get.probe <- function(name) {
  .qc.test()
  probes <- get("probes",.qcEnv)
  probes[[name]]
}

qc.add.probe <- function(name,probeset) {
  .qc.set.empty(FALSE)
  probes <- qc.get.probes()
  names  <- names(probes)
  if(name %in% names) {
    probes[name] <- probeset
  }
  else {
    probes <- c(probes,probeset)
    names  <- c(names,name)
    names(probes) <- names
  }
  assign("probes",probes,envir=.qcEnv)
}

qc.get.ratios <- function() {
  get("ratios",.qcEnv)
}

qc.get.ratio <- function(name) {
  .qc.test()
  probes <- get("ratios",.qcEnv)
  probes[[name]]
}

qc.add.ratio <- function(name,probeset1,probeset2) {
  .qc.set.empty(FALSE)
  ratios <- qc.get.ratios()
  names  <- names(ratios)
  if(name %in% names) {
    ratios[[name]] <- c(probeset1,probeset2)
  }
  else {
    if(length(ratios) == 0) {
      ratios <- list(c(probeset1,probeset2))
    }
    else {
      ratios <- c(ratios,list(c(probeset1,probeset2)))
    }
    names  <- c(names,name)
    names(ratios) <- names
  }
  assign("ratios",ratios,envir=.qcEnv)
}

qc.ok <- function() {
  !(.qc.is.empty())
}

.qc.test <- function() {
  if(.qc.is.empty()) {
    stop("QC environment is empty.\n See 'setQCEnvironment' for more details.\n")
  }
}

# read a file and use this to set up the QC Environment
# Expects a (whitespace delimited) file set out like this:
#
#
#    array <arrayname>
#    ratio <rationame> <pset1> <pset2>
#    ratio <rationame> <pset1> <pset2>
#    ratio <rationame> <pset1> <pset2>
#    ...
#    spk <spikename> <pset>
#    spk <spikename> <pset>
#    spk <spikename> <pset>
#    spk <spikename> <pset>
#    ...
#    alpha1 <value>
#    alpha2 <value>
#
# where spikename can be one of bioB, bioC, bioD or creX

qc.read.file <- function(fn) {
  .createEmptyQCEnvironment()
  fl <- file(fn,"r")
  if(!isOpen(fl)) { stop(paste("Couldn't open file",fn,"\n.")) }
  lines <- readLines(fl)
  lines <- strsplit(lines,"\\s+")
  for(l in lines) {
    a <- NULL
    if(length(l) > 0) {
      switch(l[1],
             array  = .qc.file.setQCArray(l),
             spk    = .qc.file.addQCSpike(l,a),
             ratio  = .qc.file.addQCRatio(l,a),
             alpha1 = .qc.file.setAlpha1(l,a),
             alpha2 = .qc.file.setAlpha2(l,a))
    }
  }
  close(fl)
}

.qc.file.setQCArray <- function(toks) {
  if(length(toks) != 2) {
    stop("Array name should be a single string.");
  }
  qc.set.array(toks[2])
}

#Adds or replaces the probeset for the specified QC spike on the array, 'a'
.qc.file.addQCSpike <- function(toks,a) {
  if(length(toks) != 3) {
    stop(paste("Error parsing file\nExpecting: spk '[bioB|bioC|bioD|creX] <probesetid>'\nGot: '",paste(toks,collapse=" "),"'.\n"));
  }
  if((toks[2] != "bioB") &
     (toks[2] != "bioC") &
     (toks[2] != "bioD") &
     (toks[2] != "creX"))  {
    stop(paste("Error parsing file\nSpike name must be one of  'bioB', 'bioC', 'bioD' or 'creX'\nGot: '",paste(toks,collapse=" "),"'.\n"));
  }
  qc.add.spike(toks[2],toks[3])
}

#store in the ratios list in .qcEnv each pairwise comparison. Add the relevant probesets to qc.probes, if they are not already there
.qc.file.addQCRatio <- function(toks,a) {
  if(length(toks) != 4) {
    stop(paste("Error parsing file\nExpecting: ratio <probeset1name>/<probeset2name> <probesetid1> <probesetid2>'\nGot: ",c(toks),"'.\n"));
  }
  name <- toks[2]
  ps1 <- toks[3]
  ps2 <- toks[4]
  probenames <- strsplit(name,"/")[[1]]
  if(length(probenames)!=2) {
    stop(paste("Error parsing file\nExpecting rationame: '<probeset1name>/<probeset2name>'\nGot: ",name,"'.\n"));
  }
  qc.add.probe(probenames[1],ps1)
  qc.add.probe(probenames[2],ps2)
  qc.add.ratio(name,ps1,ps2)
}

.qc.file.setAlpha1 <- function(toks,a) {
  if(length(toks) != 2) {
    stop(paste("Error parsing file\nExpecting: alpha1 <value>'\nGot: ",paste(toks),"'.\n"));
  }
  v <- as.numeric(toks[2])
  if(is.na(v)) {
    stop(paste("Error parsing file\nalpha1 value must be a number'\nGot: ",paste(toks),"'.\n"));
  }
  
  qc.set.alpha1(v)
}


.qc.file.setAlpha2 <- function(toks,a) {
  if(length(toks) != 2) {
    stop(paste("Error parsing file\nExpecting: alpha2 <value>'\nGot: ",paste(toks),"'.\n"));
  }
  v <- as.numeric(toks[2])
  if(is.na(v)) {
    stop(paste("Error parsing file\nalpha1 value must be a number'\nGot: ",paste(toks),"'.\n"));
  }
  
  qc.set.alpha2(v)
}

.qcEnv <- new.env(parent=emptyenv())

#initialize the QC environment
.createEmptyQCEnvironment <- function() {
  assign("empty",TRUE,envir=.qcEnv)
  assign("array",NULL,envir=.qcEnv)
  assign("alpha1",NULL,envir=.qcEnv)
  assign("alpha2",NULL,envir=.qcEnv)
  assign("spikes",vector(),envir=.qcEnv)
  assign("probes",vector(),envir=.qcEnv)
  assign("ratios",list(),envir=.qcEnv)
}


qc.have.params <- function(name) {
  fn <- .get.array.def.file(name)
  !(fn == "")
}


.initializeQCEnvironment <- function() {
  .createEmptyQCEnvironment();
}

.initializeQCEnvironment();


# set up the QC environment for a given array type
# if path specified, look there. Otherwise
# look in extdata to see if there's a file there
setQCEnvironment <- function(array,path=NULL) {
  if(is.null(path)) {
    fn <- .get.array.def.file(array)
    if(fn != "") {
      qc.read.file(fn)
    }
    else {
      stop(paste("Could not find array definition file '",paste(array,"qcdef",sep="."),"'. Simpleaffy does not know the QC parameters for this array type.\nSee the package vignette for details about how to specify QC parameters manually.\n"))
    }
  }
  else {
    fn <- file.path(path,paste(array,"qcdef",sep="."))
    if(fn != "") {
      qc.read.file(fn)
    }
    else {
      stop(paste("Could not find array definition file: '",paste(array,"qcdef",sep="."),"' in the directory specified:",path,"\n"))
    }
  }
}

.get.array.def.file <- function(array) {
  fn <- paste(array,"qcdef",sep=".")
  system.file("extdata",fn,package="simpleaffy")
}



.getRatios <- function(x) {
  vals <- qcProbes(x);
  if(!qc.ok()) {
    setQCEnvironment(arrayType(x))
  }
  to.calculate <- qc.get.ratios();
  r <- sapply(to.calculate,function(a) { vals[,a[1]] - vals[,a[2]]})
  return(r)
}

setGeneric("ratios", function(object) standardGeneric("ratios"))
setMethod("ratios","QCStats",function(object) .getRatios(object))


qc.affy <- function(unnormalised,normalised=NULL,tau=0.015,logged=TRUE,cdfn=cdfName(unnormalised)) {
  verbose <- getOption("verbose")
  if(.qc.is.empty()) {
    cdfn <- cleancdfname(cdfn)
    
    
    if(verbose){cat(paste("Looking cdf file for:",cdfn,"\n"))}
    
    setQCEnvironment(cdfn)
  }
  else {
    if(qc.get.array() !=  cleancdfname(cdfn)) {
      warning(paste("CDF Environment name '", qc.get.array(), "' does not match cdfname '", cleancdfname(cdfn),"'"))
    }
  }
  if(is.null(normalised)) {
    if(verbose){cat(paste("Preprocessing expression data using mas5\n"))}
    normalised <- call.exprs(unnormalised,"mas5");
  }
  
  x <- exprs(normalised);
  
  det <- detection.p.val(unnormalised,tau=tau,alpha1=qc.get.alpha1(),alpha2=qc.get.alpha2());
  
  #change
  dpv<-apply(det$call,2,function(x) {
    x[x!="P"] <- 0;
    x[x=="P"] <- 1;
    x<-as.numeric(x);
    return(100 * sum(x)/length(x));
  });
  
  sfs    <- experimentData(normalised)@preprocessing$sfs;
  target <- experimentData(normalised)@preprocessing$tgt;
  
  if(!logged) { x <- log2(x); }
  
  bgsts<-.bg.stats(unnormalised)$zonebg
  
  meanbg <- apply(bgsts,1,mean);
  
  minbg  <- apply(bgsts,1,min);
  
  maxbg  <- apply(bgsts,1,max);
  
  stdvbg <- sqrt(apply(bgsts,1,var));
  
  
  
  
  #get the probenames for the QC probes for this chip
  qc.probenames <- qc.get.probes();
  
  qc.probe.vals <- rbind(c(),(sapply(qc.probenames, function(y) {x[y,]})))
  rownames(qc.probe.vals) <- colnames(x);
  colnames(qc.probe.vals) <- qc.probenames
  spike.probenames <- qc.get.spikes();
  spike.vals <- rbind(c(),(sapply(spike.probenames, function(y) {x[y,]})))
  
  rownames(spike.vals) <- colnames(x);
  colnames(spike.vals) <- spike.probenames
  
  bb <- spike.probenames["bioB"]
  
  if(!is.na(bb)) {
    biobcalls <- det$call[bb,]
  }
  else {
    biobcalls <- NULL
  }
  
  return(new("QCStats",scale.factors=sfs,target=target,percent.present=dpv,average.background=meanbg,minimum.background=minbg,maximum.background=maxbg,
             spikes=spike.vals,qc.probes=qc.probe.vals,bioBCalls=biobcalls,arraytype=cdfn));
}


setGeneric("qc", function(unnormalised,...) standardGeneric("qc"))
setMethod("qc","AffyBatch",function(unnormalised,...) qc.affy(unnormalised,...))



.bg.stats <- function(unnormalised, grid=c(4,4)) {
  pms         <- unlist(pmindex(unnormalised))
  mms         <- unlist(mmindex(unnormalised))
  all         <- c(pms,mms)
  intensities <- exprs(unnormalised)
  rws <- nrow(unnormalised)
  cls <- ncol(unnormalised)
  zonebg <- c();
  zonesd <- c();
  for(no in 1:length(unnormalised)){
    this.array <- intensities[,no];
    result <- .C("bgmas",as.integer(as.vector(all)),as.integer(length(all)),
                 as.double(as.vector(this.array)),as.integer(length(this.array)),
                 as.integer(rws),
                 as.integer(cls),
                 as.integer(grid[1]),as.integer(grid[2]),
                 zonebg=double(grid[1] * grid[2]),
                 zonesd=double(grid[1] * grid[2]),corrected=double(length(this.array)),PACKAGE="simpleaffy");
    zonesd <- rbind(zonesd, result$zonesd);
    zonebg <- rbind(zonebg, result$zonebg);
  }
  colnames(zonesd) <- paste("zone",1:16,sep=".");
  colnames(zonebg) <- paste("zone",1:16,sep=".");
  rownames(zonesd) <- sampleNames(unnormalised);
  rownames(zonebg) <- sampleNames(unnormalised);
  return(list(zonebg=zonebg,zonesd=zonesd))
}


plot.qc.stats<-function(x,fc.line.col="black",sf.ok.region="light blue",chip.label.col="black",sf.thresh = 3.0,gdh.thresh = 1.25,ba.thresh = 3.0,present.thresh=10,bg.thresh=20,label=NULL,main="QC Stats",usemid=F,spread=c(-8,8),cex=1,...) {
  old.par <- par()
  par(mai=c(0,0,0,0))
  sfs    <- log2(sfs(x))
  
  n      <- length(sfs)
  
  meansf <- mean(sfs)
  
  dpv <- percent.present(x)
  dpv <- (round(100*dpv))/100;
  
  abg <- avbg(x)
  abg <- (round(100*abg))/100;
  
  if(is.null(label)) { label <- names(maxbg(x)) }
  d1 <- 0.0;
  d2 <- 0.0;
  d3 <- 0.0;
  
  for(i in 1:n) {
    for(j in 1:n) {
      d1 <- max(abs(sfs[i] - sfs[j]),d1);
      d2 <- max(abs(dpv[i] - dpv[j]),d2);
      d3 <- max(abs(abg[i] - abg[j]),d3);
    }
  }
  
  # set up plotting area - a column for array names next to a column for the QC
  
  m <- matrix(c(4,2,1,3) ,nrow=2,ncol=2)
  layout(m,c(1,2),c(0.1,1))
  # the title
  if(is.null(main)) { main="" }
  plot(0,0,xlim=range(0,1),ylim=range(0,1),type="n",yaxs="i",xaxt="n",yaxt="n",bty="n")
  text(0.5,0.5,labels=main,adj=0,cex=cex*2)
  
  # write out the array names
  
  plot(0,0,xlim=range(0,1),ylim=range(-1,n),type="n",yaxs="i",xaxt="n",yaxt="n",bty="n")
  text(1,(1:n)-0.5,labels=label,adj=1,cex=cex)
  plot(0,0,xlim=spread,ylim=c(-1,n),type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",bty="n")
  
  x1 <- (sf.thresh/2.0 +  meansf)
  y1 <- 0
  x2 <- (-sf.thresh/2.0 +  meansf)
  y2 <- n
  
  polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col=sf.ok.region,border=sf.ok.region);
  lines(c(0,0),c(0,n),lty=1,col=fc.line.col)
  lines(c(-1,-1),c(0,n),lty=2,col="grey")
  lines(c(-2,-2),c(0,n),lty=2,col="grey")
  lines(c(-3,-3),c(0,n),lty=2,col=fc.line.col)
  lines(c(1,1),c(0,n),lty=2,col="grey")
  lines(c(2,2),c(0,n),lty=2,col="grey")
  lines(c(3,3),c(0,n),lty=2,col=fc.line.col)
  text(3,-1,"3",pos=3,col=fc.line.col,cex=cex)
  text(2,-1,"2",pos=3,col=fc.line.col,cex=cex)
  text(1,-1,"1",pos=3,col=fc.line.col,cex=cex)
  text(-3,-1,"-3",pos=3,col=fc.line.col,cex=cex)
  text(-2,-1,"-2",pos=3,col=fc.line.col,cex=cex)
  text(-1,-1,"-1",pos=3,col=fc.line.col,cex=cex)
  text(0,-1,"0",pos=3,col=fc.line.col,cex=cex)
  
  rats <- ratios(x);
  if(!usemid) {
    gdh <- rats[,3];
    ba  <- rats[,1];
  }
  else {
    gdh <- rats[,4];
    ba  <- rats[,2];
  }
  
  bb  <- x@bioBCalls
  
  for(i in 1:n) {
    x1<-spread[1]
    x2<-spread[2]
    y1<-i-1;
    y2<-i-1;
    lines(c(x1,x2),c(y1,y2),lty=2,col="light grey")
    if(d1 > sf.thresh) { col = "red" } else {col="blue"}
    x1 <- sfs[i]
    y1 <- i-0.25
    lines(c(0,x1),c(y1,y1),col=col);
    
    points(x1,y1,col=col,pch=20);
    x2 <- gdh[i]
    y2 <- i-0.5;
    if(gdh[i] > gdh.thresh) { col = "red" } else {col="blue"}
    points(x2,y2,pch=1,col=col);
    
    x2 <- ba[i];
    y2 <- i-0.5;
    if(ba[i] > ba.thresh) { col = "red" } else {col="blue"}
    points(x2,y2,pch=2,col=col);
    
    if(d2 > present.thresh) { col = "red" } else {col="blue"}
    x2 <- spread[1]
    y2 <- i-0.25
    dpvs<-paste(dpv[i],"%",sep="")
    text(x2,y2,label=dpvs,col=col,pos=4,cex=cex);
    if(d3 > bg.thresh) { col = "red" } else {col="blue"}
    x2 <- spread[1]
    y2 <- i-0.75
    text(x2,y2,label=abg[i],col=col,pos=4,cex=cex);
    if(bb[i]!="P") {
      text(0,i-1,label="bioB",col="red",cex=cex);
    }
    
  }
  plot(0,0,xlim=range(0,1),ylim=range(0,1),type="n",yaxs="i",xaxt="n",yaxt="n",bty="n")
  if(!usemid) {
    points(0.25,0.25,pch=1)
    text(0.3,0.25,colnames(rats)[3],pos=4,cex=cex)
    points(0.25,0.5,pch=2)
    text(0.3,0.5,colnames(rats)[1],pos=4,cex=cex)
  }
  else {
    points(0.25,0.25,pch=1)
    text(0.3,0.25,colnames(rats)[4],pos=4,cex=cex)
    points(0.25,0.5,pch=2)
    text(0.3,0.5,colnames(rats)[2],pos=4,cex=cex)
  }
  
  ow <- options("warn")$warn
  options(warn=-1)
  par(old.par)
  options(warn=ow)
}

setGeneric("plot")

setMethod("plot",c("QCStats","missing"),function(x,y,...) plot.qc.stats(x,...))



#old, DEPRECATED, WILL GO SOON!

getTao <- function(name) {
  .Deprecated("qc.get.tao","simpleaffy","This function (getTao) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name);
  qc.get.tau()
}

getAlpha1 <- function(name) {
  .Deprecated("qc.get.alpha1","simpleaffy","This function (getAlpha1) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name);
  qc.get.alpha1()
}

getAlpha2 <- function(name) {
  .Deprecated("qc.get.alpha2","simpleaffy","This function (getAlpha2) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name);
  qc.get.alpha2()
}

getActin3 <- function(name) {
  .Deprecated("qc.get.ratios","simpleaffy","This function (getActin3) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name);
  rats <- qc.get.ratios()
  toks <- unlist(strsplit(names(rats),"/"))
  rats <- unlist(rats)
  j <- cbind(toks,rats)
  r <- (j[j[,1] == "actin3",2])[1]
  names(r) <- NULL
  r
}

getActinM <- function(name) {
  .Deprecated("qc.get.ratios","simpleaffy","This function (getActinM) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name);
  rats <- qc.get.ratios()
  toks <- unlist(strsplit(names(rats),"/"))
  rats <- unlist(rats)
  j <- cbind(toks,rats)
  r <- (j[j[,1] == "actinM",2])[1]
  names(r) <- NULL
  r
}

getActin5 <- function(name) {
  .Deprecated("qc.get.ratios","simpleaffy","This function (getActin5) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name);
  rats <- qc.get.ratios()
  toks <- unlist(strsplit(names(rats),"/"))
  rats <- unlist(rats)
  j <- cbind(toks,rats)
  r <- (j[j[,1] == "actin5",2])[1]
  names(r) <- NULL
  r
}

getGapdh3 <- function(name) {
  .Deprecated("qc.get.ratios","simpleaffy","This function (getGapdh3) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name);
  rats <- qc.get.ratios()
  toks <- unlist(strsplit(names(rats),"/"))
  rats <- unlist(rats)
  j <- cbind(toks,rats)
  r <- (j[j[,1] == "gapdh3",2])[1]
  names(r) <- NULL
  r
}

getGapdhM <- function(name) {
  .Deprecated("qc.get.ratios","simpleaffy","This function (getGapdhM) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name);
  rats <- qc.get.ratios()
  toks <- unlist(strsplit(names(rats),"/"))
  rats <- unlist(rats)
  j <- cbind(toks,rats)
  r <- (j[j[,1] == "gapdhM",2])[1]
  names(r) <- NULL
  r
}

getGapdh5 <- function(name) {
  .Deprecated("qc.get.ratios","simpleaffy","This function (getGapdh5) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name);
  rats <- qc.get.ratios()
  toks <- unlist(strsplit(names(rats),"/"))
  rats <- unlist(rats)
  j <- cbind(toks,rats)
  r <- (j[j[,1] == "gapdh5",2])[1]
  names(r) <- NULL
  r
}

getAllQCProbes <- function(name) {
  .Deprecated("qc.get.probes","simpleaffy","This function (getAllQCProbes) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name)
  qc.get.probes()
}

getBioB <- function(name) {
  .Deprecated("qc.get.spikes","simpleaffy","This function (getBioB) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name)
  qc.get.spikes()["bioB"]
}


getBioC <- function(name) {
  .Deprecated("qc.get.spikes","simpleaffy","This function (getBioC) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name)
  qc.get.spikes()["bioC"]
}


getBioD <- function(name) {
  .Deprecated("qc.get.spikes","simpleaffy","This function (getBioD) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name)
  qc.get.spikes()["bioD"]
}


getCreX <- function(name) {
  .Deprecated("qc.get.spikes","simpleaffy","This function (getCreX) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name)
  qc.get.spikes()["creX"]
}


getAllSpikeProbes <- function(name) {
  .Deprecated("qc.get.spikes","simpleaffy","This function (getAllSpikeProbes) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name)
  qc.get.spikes()
}


haveQCParams <- function(name) {
  .Deprecated("qc.get.spikes","simpleaffy","This function (haveQCParams) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  qc.have.params(name)
}
