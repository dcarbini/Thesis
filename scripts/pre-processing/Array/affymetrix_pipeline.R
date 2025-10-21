#Load libraries
library(zeallot)
library(tidyr)
library(readxl)
library(dplyr)
library(tools)
library(R.utils)
library(Biobase)
library(limma)
library(sva)
library(swamp)
library(infotheo)
library(biomaRt)
library(readr)
library(AnnotationDbi)
library(affy)
library(affyio)


verbose = TRUE
gse="GSE17624"


working_directory = "/nasdata/dcarbini/"
source(paste(working_directory,gse,"/scripts/functions.R",sep=""))

### set parameters
inputfile <- paste(working_directory,gse,"/metadata/updated_GSE17624.xlsx", sep="") #phenodata
fileNameColName <-  "supplementary_file_1"  #colname where cel filenames are listed -> filename .cel
sampleColName <- "geo_accession"  #GSM
arrType = "af_exp"
dyeColName <- NULL  #colnames where dye names are listed  riflessioni123

#-------------------------------------------------------------------------- Metadata -------------------------------------------------------------------------------------------
result <- load_phenotype(inputfile)
orig_phTable <- result[[1]] #different format
phTable <- result[[2]]
phClassTable <- result[[3]] #types of variables 
phFactor <- convert_to_chr(phTable)
#---------------------------------------------------------------------------Raw data ---------------------------------------------------------------------------------------
### LOAD RAWDATA
#params for loading rawdata
removedSamplesInfo <- c()   #row of samples to remove from analysis
celDir <- paste(working_directory, gse, "/raw_data/",sep="") #directory raw_data (no .tar)
affCDF <- "hgu133plus2cdf" #to be changed depending on the array

#path where cel files are stored
#needCDF= F #bool to identify whether affy platform is the one needing cdf files or not
affyBatchObject <- NULL

#clean the tables from useless information
clean_phTable <- as.data.frame(phTable)
orig_phTable <- as.data.frame(orig_phTable)
phFactor <- as.data.frame(phFactor)

if (arrType == "af_exp") {
  results <- upload_rawdata(clean_phTable, orig_phTable, phFactor, removedSamplesInfo, fileNameColName, sampleColName, dyeColName, celDir, arrType, affCDF)
  # Now data are already normalized  
}

norm.data = as.data.frame(results$exprs)
#columns are the samples, rows are the probes
boxplot(norm.data, main = "after normalization", las=2, cex=0.7,col = rgb(1, 0, 0, alpha = 0.5))


#----------------------------------------------------------------------- Annotation trough CDF -----------------------------------------------------------------------
# hgu133plus2.db
# Estract ENSEMBLID and gene symbol
annotat <- AnnotationDbi::select(hgu133plus2.db,               
                                 keys = rownames(norm.data),  
                                 columns = c("ENSEMBL", "SYMBOL"),  
                                 keytype = "PROBEID")


norm.data$PROBEID <- rownames(norm.data)
result = merge(norm.data, annotat, by="PROBEID")

result$PROBEID <- c()
result$SYMBOL <- c()

#Since many probes don't map to ensembl gene filter them out
norm_data_filtered <- result %>%
  filter(!is.na(ENSEMBL))

# If there are probes that map to same ENSEMBL ID take the median of the expression value
# In affymetrix you have 2 probes on the same gene, so to have 1 value you take the median
selected_cols = setdiff(colnames(norm_data_filtered), c("PROBEID","ENSEMBL","SYMBOL"))
norm_data_median <- aggregate(norm_data_filtered[,selected_cols], 
                              by = list(ENSEMBL = norm_data_filtered$ENSEMBL), 
                              FUN = function(x) median(x, na.rm = TRUE))
norm_data_median$ENSEMBL -> rownames(norm_data_median)
norm_data_median$ENSEMBL <- c()

write.table(norm_data_median,file=paste(working_directory,gse, "/results2/normalized_expression_matrix_",gse,".txt",sep=""),quote = FALSE, sep = "\t", row.names = TRUE, col=NA)

norm.data <- as.matrix(norm_data_median)
rownames_preserved <- rownames(norm.data)
norm.data <- apply(norm.data, 2, as.numeric)
rownames(norm.data) <- rownames_preserved


#BEFORE TO PROCEED WITH MODELING AND BATCH EFFECTS ASSESSMENT, DEVDE THE PD AND THE NORMALIZED MATRIX BY THE TIME
# If you have different cell_lines in the same GSE is easier to divide them in different datasets. You have to divide both the phenodata (using the 'cell_line' column as discriminant) and the normalized matrix (using the GSM on the column as discriminant). 
# If you have different time points you can divide more looking also at it. The aim of this procedure is to separate the controls in the correct groups.


# conditions to divide the dataset
conds <- c("8h", "24h", "48h")

# loop 
for (cond in conds) {
  
  # Filtro phenotype e normalizzati per la condizione
  pd <- clean_phTable[clean_phTable$exposure_time == cond, ]
  column_subset <- pd$geo_accession
  norm.data.sub <- norm.data[, column_subset]
  
  #----------------------------------------------------MODELING BASED ON BATCH CORRECTION OF UNKNOWN VARIABLES---------------------------------------------------------
  
  # preparing design based on conditions and possible surrogate variables # https://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf
  # In the pd have to be present only the used conditions, so if you have 2 kind of control you have to delete the one that you don't use 
  
  # Building the group
  pd$group <- paste(pd$exposure, pd$"dose(amount)", pd$"dose(unit)", pd$exposure_time, sep = "_")
  
  # Models
  mod <- model.matrix(~0 + group, data = pd)# full model matrix - including both the adjustment variables and the variable of interest
  mod0 <- model.matrix(~1, data = pd)# the null model contains only the adjustment variables. Since we are not adjusting for any other variables in this analysis, only an intercept is included in the model.
  
  # adding surrogate variables to design matrix
  svobj <- sva(dat = norm.data.sub, mod = mod, mod0 = mod0, n.sv = NULL)# n.sv null means number of latent factors will be estimated for you
  # with the sva you have the number of covariates 
  modSv <- cbind(mod, svobj$sv)# including variable of interest
  mod0Sv <- cbind(mod0, svobj$sv) # does not include vairable of interest
  # we use the matrices with the covariates even if we don't correct for them 
  
  # Linear model fitted using covariates
  fit <- lmFit(norm.data.sub, modSv)
  
  # adjusting conditions for further comparison in contrast
  conditions <- unique(pd$group)# get the unique conditions (levels of exposure)
  colnames(mod) <- conditions
  
  # if covariates are found, rename matrix
  if(length(conditions) < ncol(modSv)){
    n_cov = abs(length(conditions) - ncol(modSv))
    covariate_names <- paste0("cov", seq_len(n_cov))
    new_conditions <- c(conditions, covariate_names)
  }
  colnames(modSv) <- new_conditions
  
  # change the phenodata to add the covariates
  idx = grepl(pattern="cov",colnames(modSv))
  cov_model = modSv[,idx]
  pd_mod = cbind(pd, cov_model); pd_mod = as.data.frame(pd_mod)
  
  write.table(pd_mod, file = paste0(working_directory, gse, "/results2/model_cov_", gse, "_", cond, ".txt"),
              quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
  
  
  #--------------------------------------------------------------- Differential expression analysis --------------------------------------------------------------------
  
  # Generation of dynamic contrasts based on the variables of interest
  contrast_list <- c()
  for (i in 2:length(conditions)) {  
    contrast <- paste(conditions[i], "-", conditions[1],sep="")
    contrast_list <- c(contrast_list, contrast)
  }
  #conditions[1] is the control -> change if in your data is different
  
  # to make all pair-wise comparisons between the groups, we do not include the surrogate variables in the contrasts (since they are only being used to adjust the analysis)
  contrast.matrix <- makeContrasts(contrasts = contrast_list, levels = modSv) # covariates are being included as part of the model, but they are not part of the contrast.
  colnames(fit$coefficients) <- new_conditions
  fit2 <- eBayes(contrasts.fit(fit, contrast.matrix))# correct looking at the contrast list
  
  # Extraction differentially expressed genes (DEGs)
  list.top.tables <- list()
  for(i in 1:length(contrast_list)) {
    list.top.tables[[i]] <- topTable(fit2, coef=i, p.value=1, lfc=0, 
                                     adjust.method="BH",
                                     number=nrow(norm.data), sort.by="logFC")
  }#p.value=1, lfc=0 = no filtering -> take everything
  names(list.top.tables) <- contrast_list
  
  # create the dataframe
  degs_all_conditions <- data.frame()
  for (comparison in seq_along(contrast_list)) {
    result <- list.top.tables[[comparison]]
    result$EnsemblID <- rownames(result)
    rownames(result) <- NULL
    group_comparison <- contrast_list[[comparison]]
    group <- unlist(strsplit(group_comparison, split = "-"))[[1]]
    degs_all_conditions <- rbind(degs_all_conditions,
                                 cbind(result, group_comparison = group_comparison, group = group))
  }
  degs_all_conditions <- unique(degs_all_conditions)
  
  write.csv(degs_all_conditions, file = paste0(working_directory, gse, "/results2/", gse, "_", cond, "_unfiltered_DEGs.csv"))
  
  # Filtering
  #cut off lfc>=0.58 and adj.p.val<=0.05 
  degs_filtered <- degs_all_conditions[
    (degs_all_conditions$logFC >= 0.58 | degs_all_conditions$logFC <= -0.58) & 
      degs_all_conditions$adj.P.Val <= 0.05, 
    ]
  write.csv(degs_filtered, file = paste0(working_directory, gse, "/results2/", gse, "_", cond, "_filtered_DEGs.csv"))
  

  # Confounding analysis
  test <- sapply(colnames(pd_mod), function(b) length(table(pd_mod[, b])) > 1)
  pd_mod[] <- lapply(pd_mod, function(x) if (is.character(x)) as.factor(x) else x)
  
  pdf(paste0(working_directory, gse, "/results2/confounding_", cond, ".pdf"), width = 8, height = 8)
  confounding_results <- swamp::confounding(pd_mod[, test], margins = c(10, 10))
  dev.off()
  
  rownames(pd_mod) <- pd_mod$geo_accession
  pr <- swamp::prince(norm.data.sub, pd_mod[, test], top = 5)
  pdf(paste0(working_directory, gse, "/results2/prince_plot_", cond, ".pdf"), width = 15, height = 15)
  swamp::prince.plot(prince = pr, margins = c(15, 15))
  dev.off()
  
  swamp::hca.plot(norm.data.sub, pd_mod[, test], method = "correlation")
  
  
}

