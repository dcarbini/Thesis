library(purrr)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(AOPfingerprintR)
library(writexl)
library(readxl)
library(dplyr)
library(tidyr)
library(circlize) 
library(stringr)


################### Combine KE uniques ############################
folder <- "C:/Users/Utente/Desktop/denys/results/"

files <- c(
  "GSE211183" = paste0(folder, "GSE211183/results/bmdx/phenodata_bmdx_noCyt.xlsx"),
  "GSE50705" = paste0(folder, "GSE50705/results/bmdx/phenodata_bmdx.xlsx"),
  "GSE17624" = paste0(folder, "GSE17624/results2/bmdx/log/phenodata_bmdx.xlsx"),
  "GSE86472" = paste0(folder, "GSE86472/results/bmdx/phenodata_bmdx.xlsx"),
  "GSE153320" = paste0(folder, "GSE153320/results/bmdx/log_nocyt/phenodata_bmdx_noCyt.xlsx"),
  "GSE69851_HepG2" = paste0(folder, "GSE69851/HepG2/bmdx/phenodata_bmdx.xlsx"),
  "GSE69851_HepaRG" = paste0(folder, "GSE69851/HepaRG/bmdx/phenodata_bmdx.xlsx"),
  "GSE69851_MCF7" = paste0(folder, "GSE69851/MCF7/bmdx/phenodata_bmdx.xlsx"),
  "GSE249377_24h_CS" = paste0(folder, "GSE249377/results/24h_CS/not_divided/bmdx/phenodata_bmdx.xlsx"),
  "GSE249377_24h_HI" = paste0(folder, "GSE249377/results/24h_HI/not_divided/bmdx/phenodata_bmdx.xlsx"),
  "GSE249377_6h_CS" = paste0(folder, "GSE249377/results/6h_CS/not_divided/bmdx/phenodata_bmdx.xlsx"),
  "GSE249377_12h_CS" = paste0(folder, "GSE249377/results/12h_CS/not_divided/bmdx/phenodata_bmdx.xlsx"),
  "GSE249377_12h_HI" = paste0(folder, "GSE249377/results/12h_HI/not_divided/bmdx/phenodata_bmdx.xlsx")
)

find_min_max_ke <- function(f){
  sheet_names <- excel_sheets(f)
  result_df <- map_dfr(sheet_names, function(sheet) {
    df <- read_excel(f, sheet = sheet)
    
    # compute min and max
    min_val <- min(df$dose_amount[df$dose_amount != 0], na.rm = TRUE)
    max_val <- max(df$dose_amount, na.rm = TRUE)
    range_val <- max_val-min_val
    
    # put everything in a dataframe
    tibble(
      Experiment = sheet,
      min = min_val,
      max = max_val,
      range = range_val
    )
  })
}

all_minmax_ke <- do.call(rbind, Map(function(f, name) {
  df <- find_min_max_ke(f)
  df$Dataset <- name
  df
}, files, names(files)))
row.names(all_minmax_ke) <- NULL
all_minmax_ke$id <- paste(all_minmax_ke$Dataset, all_minmax_ke$Experiment, sep="_")

write_xlsx(all_minmax_ke, paste0(folder,"ke_enrichment/min_max_ke.xlsx"))


#Set the directory
#dir_path <- "C:/Users/carbi/Desktop/ke_enrichment/"
dir_path <- "C:/Users/Utente/Desktop/denys/results/ke_enrichment/"

#Function to take only the unique Events in a dataset
unique_events <- function(df) {
  df <- df %>% distinct(TermID, Experiment, .keep_all = TRUE) 
  return(df)
}

#Function to map the dataframe with the Biological_system_annotations
merge_annotation <- function(df){
  #load("C:/Users/carbi/Downloads/Biological_system_annotations.rda")
  #load("C:/Users/Utente/Desktop/denys/scripts/bmdx_noshiny-master/Biological_system_annotations.rda")
  annot <- Biological_system_annotations
  # Rename one of the 2 columns 
  names(annot)[names(annot) == "ke"] <- "TermID"
  
  # Merging on the dataframes on the column "TermId"
  final_df <- merge(df, annot, by = "TermID")
  return(final_df)
}

normalize_bmd_ke <- function(df, minmax){
  minmax <- minmax %>% select(min, max, range, id)
  
  #Force the column as numeric
  df <- df %>%
    mutate(BMD = as.numeric(BMD),
           id = paste(Dataset, chemical, sep="_")
    )%>%
    left_join(minmax, by = "id")
  
  df_normalized <- df %>%
    mutate(BMD_norm = ifelse(range == 0,
                             1.5,
                             1 + (BMD - min) / range))
  return(df_normalized)
}


#Function to read and combine excel files with the same name
#NOrmalizaition 
combine_excels <- function(pattern, path = dir_path) {
  files <- list.files(path = path, pattern = paste0("^", pattern), full.names = TRUE)
  df_list <- lapply(files, function(file) {
    df <- read_excel(file)
    df$Dataset <- str_extract(basename(file), "GSE[^.]+")
    #Call the function to take only the unique events
    df <- unique_events(df)
    #Call the function to merge the dataframe with the Biological_system_annotations
    df <- merge_annotation(df)%>% 
      separate(Experiment, into = c("chemical", "time"), sep = "_", remove = FALSE) %>%
      mutate(time = as.numeric(time))
    #Normalize the BMD values
    #all_minmax_ke <- read_excel("C:/Users/Utente/Desktop/denys/results/ke_enrichment/min_max_ke.xlsx")
    df <- normalize_bmd_ke(df, all_minmax_ke)
    return(df)
  })
  combined_df <- bind_rows(df_list)
  return(combined_df)
}


# Create the dataframes
ke_df <- combine_excels("ke_enrichment_results", paste0(dir_path, "ke_enrichment_files/"))


#Checks
View(ke_df)

write_xlsx(ke_df, paste0(dir_path,"ke.xlsx"))


################### Combine KE NON uniques ############################

combine_nounique <- function(pattern, minmax_ke, path = dir_path) {
  files <- list.files(path = path, pattern = paste0("^", pattern), full.names = TRUE)
  df_list <- lapply(files, function(file) {
    df <- read_excel(file)%>% 
      separate(Experiment, into = c("chemical", "time"), sep = "_", remove = FALSE) %>%
      mutate(time = as.numeric(time))
    df$Dataset <- str_extract(basename(file), "GSE[^.]+")
    #minmax_ke <- read_excel("C:/Users/Utente/Desktop/denys/results/ke_enrichment/min_max_ke.xlsx")
    df <- normalize_bmd_ke(df, minmax_ke)
    return(df)
  })
  combined_df <- bind_rows(df_list)
  return(combined_df)
}

ke_df_nonunique <- combine_nounique("ke_enrichment_results", all_minmax_ke, paste0(dir_path, "ke_enrichment_files/"))

write_xlsx(ke_df_nonunique, paste0(dir_path,"ke_nonunique.xlsx"))

combine_excels_no_norm <- function(pattern, path = dir_path) {
  files <- list.files(path = path, pattern = paste0("^", pattern), full.names = TRUE)
  df_list <- lapply(files, function(file) {
    df <- read_excel(file)%>% 
      separate(Experiment, into = c("chemical", "time"), sep = "_", remove = FALSE) %>%
      mutate(time = as.numeric(time))
    df$Dataset <- str_extract(basename(file), "GSE[^.]+")
    return(df)
  })
  combined_df <- bind_rows(df_list)
  return(combined_df)
}

aop_df <- combine_excels_no_norm("aop_enrichment_results", paste0(dir_path, "aop_enrichment_files/"))
write_xlsx(aop_df, paste0(dir_path,"aop.xlsx"))

#################### Combine all DDGs ########################################

folder <- "C:/Users/Utente/Desktop/denys/results/"

files_ddg <- list.files(paste0(folder,"DoseDependentGenes/DDG/"), pattern = "^GSE.*\\.xlsx$", full.names = TRUE)

process_file <- function(f) {
  df <- read_excel(f, col_names = TRUE)
  
  #Map ensembl to symbols
  df <- df %>%
    left_join(genes_human, by = c("Feature" = "ensembl_gene_id")) %>%
    mutate(Feature = ifelse(!is.na(hgnc_symbol) & hgnc_symbol != "", 
                            hgnc_symbol, 
                            Feature)) %>%
    select(-hgnc_symbol)
  
  df$Dataset <- gsub("\\.[[:alnum:]]+$", "", basename(f))
  
  return(df)
}

#Unify files
all_data <- bind_rows(lapply(files_ddg, process_file))

dim(all_data)

##################### Take min and max doses per each chemical and dataset ##########################################
folder <- "C:/Users/Utente/Desktop/denys/results/"

files <- c(
  "GSE211183" = paste0(folder, "GSE211183/results/bmdx/phenodata_bmdx_noCyt.xlsx"),
  "GSE50705" = paste0(folder, "GSE50705/results/bmdx/phenodata_bmdx.xlsx"),
  "GSE17624" = paste0(folder, "GSE17624/results2/bmdx/log/phenodata_bmdx.xlsx"),
  "GSE86472" = paste0(folder, "GSE86472/results/bmdx/phenodata_bmdx.xlsx"),
  "GSE153320" = paste0(folder, "GSE153320/results/bmdx/log_nocyt/phenodata_bmdx_noCyt.xlsx"),
  "GSE69851_HepG2" = paste0(folder, "GSE69851/HepG2/bmdx/phenodata_bmdx.xlsx"),
  "GSE69851_HepaRG" = paste0(folder, "GSE69851/HepaRG/bmdx/phenodata_bmdx.xlsx"),
  "GSE69851_MCF7" = paste0(folder, "GSE69851/MCF7/bmdx/phenodata_bmdx.xlsx"),
  "GSE249377_24H_CS" = paste0(folder, "GSE249377/results/24h_CS/not_divided/bmdx/phenodata_bmdx.xlsx"),
  "GSE249377_24H_HI" = paste0(folder, "GSE249377/results/24h_HI/not_divided/bmdx/phenodata_bmdx.xlsx"),
  "GSE249377_6H_CS" = paste0(folder, "GSE249377/results/6h_CS/not_divided/bmdx/phenodata_bmdx.xlsx"),
  "GSE249377_12H_CS" = paste0(folder, "GSE249377/results/12h_CS/not_divided/bmdx/phenodata_bmdx.xlsx"),
  "GSE249377_12H_HI" = paste0(folder, "GSE249377/results/12h_HI/not_divided/bmdx/phenodata_bmdx.xlsx")
)

find_min_max <- function(f){
  sheet_names <- excel_sheets(f)
  result_df <- map_dfr(sheet_names, function(sheet) {
    df <- read_excel(f, sheet = sheet)
    
    # compute min and max
    min_val <- min(df$dose_amount[df$dose_amount != 0], na.rm = TRUE)
    max_val <- max(df$dose_amount, na.rm = TRUE)
    range_val <- max_val-min_val
      
    # put everything in a dataframe
    tibble(
      Experiment = sheet,
      min = min_val,
      max = max_val,
      range = range_val
    )
  })
}

all_minmax <- do.call(rbind, Map(function(f, name) {
  df <- find_min_max(f)
  df$Dataset <- name
  df
}, files, names(files)))
row.names(all_minmax) <- NULL
all_minmax$id <- paste(all_minmax$Dataset, all_minmax$Experiment, sep="_")

write_xlsx(all_minmax, paste0(folder,"DoseDependentGenes/min_max.xlsx"))

###################### Normalize BMD in DDGs file ####################################

nemesis <- c("BPA", "BPF", "BPS", "DEP", "DINP", "Hexamoll DINCH", "PFOA", "PFHxS", "ZEA")
eu <- c("BPA", "BPB", "BPS", "BPAF", "paranonylphenol", "Prochloraz", "Propiconazole", "Ziram", "4-(4-propan-2-yloxyphenyl)sulfonylphenol", "dibutylphthalate", "lithiumchloride", "DEHP", "triclosan", "TBBA")
zenodo <- c('RU486','Clobetasol', 'Dihydroxyvitamin3','Dehydrorepiandrosterone',  'BPA','Genistein', 'ValproicAcid', 'Tamoxifen', 'Progesterone', 'Ketoconazole', 'Phenobarbital', 'Flutamide', 'Trenbolone', 'Griseofulvin', 'Methanol','Estradiol', 'Nicotine','Vinblastine', 'FK506/Cyclosporin A',  'diethylstilbestrol', 'Daidzein', 'Ethynylestradiol', 'paranonylphenol', 'BPS', 'BPF',  'Glyphosate', 'Roundup', 'BPAF', 'BPAP','BPB', 'BPZ', 'PPT', 'benzo[a]pyrene', 'Ethanol','TBBA','pendimethalin','stomp_aqua_pendimethalin', 'Reserpine', 'Thiram', 'Propiconazole','Cyproterone.acetate', 'Bifenthrin', 'Trifloxystrobin', 'PFOA', 'Cycloheximide','Simazine', 'Cyproconazole', 'Simvastatin', 'Pyraclostrobin', 'PFOS', 'Triiodothyronine', 'Vinclozolin', '4.Hydroxytamoxifen', 'Maneb','Tetrac', 'Ziram','4.Cumylphenol', 'Lovastatin', 'Cyanazine', 'Fulvestrant', 'Nilutamide','Cypermethrin', 'Prochloraz', 'gesaprim', 'Imazalil', 'Rotenone','Dexamethasone', 'ZEA','TGSA', 'BADGE')
more <- c("2,4.BPF","2,4.BPS", "17B.estradiol", "4,4.BPF", "BPS.MPE", "Dex", "para.nonylphenol")
accepted_chemicals <- unique(c(nemesis, eu, zenodo, more))

normalize_bmd <- function(df, minmax){
  minmax$id <- paste(minmax$Dataset, minmax$Experiment, sep="_")
  minmax <- minmax %>% select(min, max, range, id)
  
  #Force the column as numeric
  df <- df %>%
    mutate(BMD = as.numeric(BMD),
           id = paste(Dataset, Experiment, sep="_")
           )%>%
    left_join(minmax, by = "id")
    
  df_normalized <- df %>%
    mutate(BMD_norm = ifelse(range == 0,
                             1.5,
                             1 + (BMD - min) / range))
  return(df_normalized)
}


all_norm <- normalize_bmd(all_data, all_minmax)%>%
  mutate(Experiment = ifelse(Experiment == "paranonylphenol", "para.nonylphenol", Experiment),
         Experiment = ifelse(Experiment == "17B.estradiol", "Estradiol", Experiment))%>%
  filter(Experiment %in% accepted_chemicals)

View(all_norm)
write_xlsx(all_norm, paste0(folder,"DoseDependentGenes/all_doseDependent.xlsx"))

################## Make a unique DEGs file ####################################


library(dplyr)
library(readr)
library(stringr)
library(purrr)

# cartella con i CSV
cartella <- "C:/Users/Utente/Desktop/denys/results/comb/DEGs/"

# leggi i file csv
file_list <- list.files(cartella, pattern = "\\.csv$", full.names = TRUE)

# funzione per processare ogni file
process_file <- function(file) {
  df <- read_csv(file, show_col_types = FALSE)%>%
    select(-starts_with("...")) %>%
    mutate(across(everything(), as.character))
  
  # se esiste la colonna EnsemblID
  if ("EnsemblID" %in% colnames(df)) {
    df <- df %>%
      left_join(genes_human, by = c("EnsemblID" = "ensembl_gene_id")) %>%
      mutate(Feature = ifelse(!is.na(hgnc_symbol) & hgnc_symbol != "",
                              hgnc_symbol,
                              EnsemblID)) %>%
      select(-hgnc_symbol)
  } else {
    colnames(df)[colnames(df) == "GeneSymbol"] <- "Feature"
  }
  df <- df %>% mutate(Feature = as.character(Feature))
  
  # aggiungi chemical da "group"
  if ("group" %in% colnames(df)) {
    df <- df %>%
      mutate(Experiment = str_split_fixed(group, "_", 2)[,1])
  } else {
    df <- df %>% mutate(Experiment = NA)
  }
  
  # aggiungi dataset dal nome file
  dataset_name <- str_replace(basename(file), "_filtered.*", "")
  df <- df %>% mutate(Dataset = dataset_name)
  
  return(df)
}

# applica a tutti i file
deg <- map_dfr(file_list, process_file)

deg <- deg %>%
  mutate(Dataset = ifelse(str_detect(Dataset, "GSE69851_HepaRG"),
                          "GSE69851_HepaRG",
                          Dataset),
         Dataset = ifelse(str_detect(Dataset, "GSE69851_MCF7"),
                          "GSE69851_MCF7",
                          Dataset),
         Dataset = ifelse(str_detect(Dataset, "GSE50705_"),
                          "GSE50705",
                          Dataset))
dict <- c("GSE153320" = "72", "GSE17624_8h" = "8", "GSE17624_24h" = "24", 
          "GSE17624_48h" = "48", "GSE211183" = "48", "GSE249377_12h_CS" = "12", 
          "GSE50705" = "Type3", "GSE249377_24h_CS" = "24", "GSE249377_12h_HI" = "12", 
          "GSE86472" = "48", "GSE69851_MCF7" = "6", "GSE69851_HepG2" = "6",
          "GSE69851_HepaRG" = "6", "GSE249377_24h_HI" = "24", "GSE249377_6h_CS" = "6")

deg <- deg %>%
  mutate(exposure_time = recode(Dataset, !!!dict))


write_xlsx(deg, paste0(folder,"DoseDependentGenes/all_DEGs.xlsx"))

##############################################################################
