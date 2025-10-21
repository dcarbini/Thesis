library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(purrr)
library(RColorBrewer)
library(factoextra)
library(writexl)

# === 1. Funzione per leggere un file Excel ===
read_aop_file <- function(file_path) {
  df <- read_excel(file_path)
  
  # Rimuovi duplicati mantenendo tutte le colonne
  df <- df %>%
    distinct(TermID, Experiment, .keep_all = TRUE)
  
  # Estrai GSE dal nome file
  gse <- str_extract(basename(file_path), "GSE[^.]+")
  df$GSE <- gse
  
  return(df)
}

read_aop_file_withduplicates <- function(file_path) {
  df <- read_excel(file_path)
  
  # Estrai GSE dal nome file
  gse <- str_extract(basename(file_path), "GSE[^.]+")
  df$GSE <- gse
  
  return(df)
}

# === 2. Caricamento di tutti i file Excel ===
# Sostituisci con il percorso della tua cartella
#folder_path <- "C:/Users/carbi/Desktop/aop_enrichment/"
folder_path <- "C:/Users/Utente/Desktop/denys/aop_enrichment/"
excel_files <- list.files(folder_path, pattern = "^aop.*\\.xlsx$", full.names = TRUE)

# Leggi tutti i file e combinali
all_data <- map_dfr(excel_files, read_aop_file)
all_duplicates <- map_dfr(excel_files,read_aop_file_withduplicates)

# === 3. Analisi ===

# Numero di dataset
n_datasets <- length(unique(all_data$GSE))

# A. Pathway condivisi
shared_pathways <- all_data %>%
  group_by(TermID, a.name) %>%
  summarise(n_datasets = n_distinct(GSE), .groups = 'drop') %>%
  filter(n_datasets > 1)

# B. Pathway unici
unique_pathways <- all_data %>%
  group_by(GSE, TermID, a.name) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(TermID) %>%
  filter(n_distinct(GSE) == 1)

# D. Matrice per clustering
# Creiamo una matrice GSE x TermID (1 se presente, 0 altrimenti)
binary_matrix <- all_data %>%
  select(GSE, TermID) %>%
  distinct() %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = TermID, values_from = value, values_fill = 0)

# Heatmap dei pathway
heatmap_data <- as.matrix(binary_matrix[,-1])
rownames(heatmap_data) <- binary_matrix$GSE
png(paste0(folder_path,"heatmap_aop.png"), width = 2500, height = 600, res = 150)
# ...Heatmap()...
Heatmap(heatmap_data, name = "Presence", col = c("white", "blue"),
        show_row_names = TRUE, show_column_names = TRUE,
        cluster_rows = TRUE, cluster_columns = TRUE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        column_names_gp = gpar(fontsize = 6),
        column_title = "Presence of AOPs in AOP fingerprints")
dev.off()

# E. Numero di pathway per GSE
pathways_per_gse <- all_data %>%
  group_by(GSE) %>%
  summarise(n_pathways = n_distinct(TermID))

ggplot(pathways_per_gse, aes(x = GSE, y = n_pathways)) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  labs(title = "Numero di Pathway per Dataset", y = "N° Pathway", x = "GSE")

# G. Esportazione risultati
write_xlsx(shared_pathways, paste0(folder_path,"shared_pathways.xlsx"))
write.csv(unique_pathways, paste0(folder_path,"unique_pathways.csv"), row.names = FALSE)



########################################

#Quali experiments mi danno un termID
term_id_target <- "Aop:393"

all_data %>%
  filter(TermID == term_id_target) %>%
  distinct(Experiment) %>%
  pull(Experiment)



