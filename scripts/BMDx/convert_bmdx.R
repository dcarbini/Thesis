# upload the libraries
library(readxl)
library(writexl)
library(dplyr)
library(stringr)

GSE <- "GSE153320"
dir <- paste0("/project/omics/public_data/endocrine_disruptors/",GSE,"/results/")

# upload the original files
meta_data <- read_excel(paste0(dir,"updated_GSE153320_DEG.xlsx"))
expr_matrix <- tryCatch({
  read_excel(paste0(dir,"vst_expression_matrix_GSE153320_Blind_F.txt"))
}, error = function(e) {
  message("NOTE: It is not an excel file, trying to read it as table...")
  read.table(paste0(dir,"vst_expression_matrix_GSE153320_Blind_F.txt"), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
})

output_dir <- paste0(dir, "bmdx/")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("\n==== Data Uploaded ====\n")

# extract the list of chemicals 
chemicals <- unique(str_split_fixed(meta_data$DEG_variables, "_", 2)[, 1])
# remove controls from the list
chemicals <- chemicals[!is.na(chemicals) & chemicals != "NA" & chemicals != "" & chemicals != "DMSO"]  

cat("\n==== Chemical extracted ====\n")
cat(paste0("In the dataset are present ", length(chemicals), " chemicals."))

# declare the outputs
phenodata_list <- list()
vst_counts_list <- list()

# iterate over the chemicals to split original dataframes
for (chem in chemicals) {
  
  cat(paste0("\n==== Processing: ",chem, " ====\n"))
  
  # subset meta_data: chemical + controls
  subset_meta <- meta_data %>%
    filter(str_detect(DEG_variables, paste0("^", chem, "_")) | str_detect(DEG_variables, "NA_") | str_detect(DEG_variables, "DMSO_"))
  
  # add to the phenodata output the subset
  phenodata_list[[chem]] <- subset_meta
  
  # subset the expression matrix: take the samples of the metadata subset looking at the geo_accession
  vst_counts_list[[chem]] <- cbind(Gene = rownames(expr_matrix), expr_matrix[, subset_meta$geo_accession, drop = FALSE])
}

# save the results
write_xlsx(phenodata_list, paste0(output_dir, "phenodata_bmdx.xlsx"))
write_xlsx(vst_counts_list, paste0(output_dir, "vst_counts_bmdx.xlsx"))